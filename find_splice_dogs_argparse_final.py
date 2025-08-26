#how to run python this.script.py project/path --run_all (can also do --run_splice to only find splice dogs)
#folder setup:
    # bed folder 
import time
import argparse
import pandas as pd 
import numpy as np
from pathlib import Path
import glob
import subprocess
import os
project_dir = Path(__file__).resolve().parent
from dogtools import extract_Pol2_genes as pol2
from dogtools.gtf import gtf_to_bed

def generate_and_validate_sample_matrix(bed_files, sample_matrix_path):
    prefixes = sorted(set([Path(f).stem.split("_")[0].split(".")[0] for f in bed_files]))
    if not os.path.exists(sample_matrix_path):
        print("Sample matrix not found. Generating template...")
        with open(sample_matrix_path, "w") as f:
            f.write("sample\tcondition\n")
            for p in prefixes:
                f.write(f"{p}\tcontrol_or_test\n")
        raise FileNotFoundError(
            f"Sample matrix created at {sample_matrix_path}. "
            f"Please fill in the 'condition' column (e.g., 'control' or 'zta') and rerun the script.")

    df = pd.read_csv(sample_matrix_path, sep="\t")
    missing = [p for p in prefixes if p not in df["sample"].values]
    if missing:
        raise ValueError(f"Sample matrix is missing entries for: {missing}")
    return sample_matrix_path

def find_position(arr, val): 
    '''Finds the position of a junction in an array of coordinates. Must be a sorted array'''

    lower_bound = 0 
    upper_bound = len(arr) 
    while lower_bound < upper_bound: 
        middle = (lower_bound + upper_bound) // 2 #defines middle section; '//2' = rounds down to nearest whole number
        if val >= arr[middle]: #looks for splice junction in middle section of annotation file, if splice junction is not within that section, update boundaries and continue search
            lower_bound = middle + 1 
        else:
            upper_bound = middle
    return lower_bound - 1 #if the lower_bound is no longer < upper_bound, then go start over with new boundaries

def summarize_splice_dogs(records, file_label):
    rows = []
    for row in records:
        gene_info = f"{row['gene_chrom']}:{row['gene_start']}-{row['gene_end']}:{row['gene_strand']}:{row['gene_name']}"
        sj_coord = f"{row['sj_chrom']}:{row['sj_start']}-{row['sj_end']}:{row['sj_strand']}"
        sj_count = row['sj_count']
        rows.append({
            "gene info": gene_info,
            "Dog_sj_coordinate": sj_coord,
            "# of sj_count from bed": sj_count,
            "file_name": file_label
        })
    return pd.DataFrame(rows)


def build_splice_dogs(bed_files, genes_only, output_dir):
    all_summaries = []
    for bed_file in bed_files:
        print(f"Processing: {bed_file}")
        file_path = Path(bed_file)
        stem = file_path.stem
        sample_part = stem.split("_")[0]
        strand_part = "+" if "positive" in stem else "-"
        file_label = f"{sample_part}:{strand_part}"

        junctions = pd.read_csv(file_path, sep='\t', header=None, comment="#")
        dog_results = []

        for chrom in genes_only["chromosome"].unique():
            if strand_part == "+":
                filtered_genes = genes_only[(genes_only["chromosome"] == chrom) & (genes_only["strand"] == '+')].sort_values('left')
            else:
                filtered_genes = genes_only[(genes_only["chromosome"] == chrom) & (genes_only["strand"] == '-')].sort_values('right')

            filtered_junctions = junctions[junctions[0] == chrom]
            if filtered_genes.empty or filtered_junctions.empty:
                continue

            gene_coords = list(filtered_genes['left'] if strand_part == "+" else filtered_genes['right'])
            filtered_genes.index = range(len(filtered_genes.index))

            for sj in filtered_junctions.index:
                sj_start = filtered_junctions.loc[sj, 1]
                sj_end = filtered_junctions.loc[sj, 2]

                if strand_part == "+":
                    index_start = find_position(gene_coords, sj_start)
                    overlap_genes = []
                    for gene in filtered_genes.index[index_start:]:
                        gene_start = filtered_genes.loc[gene, "left"]
                        gene_end = filtered_genes.loc[gene, "right"]
                        if gene_start <= sj_start <= gene_end:
                            overlap_genes.append([gene, gene_start, gene_end])
                        elif gene_start > sj_end:
                            break

                    # Find longest downstream gene
                    longest_gene = None
                    max_end = 0
                    for g in overlap_genes:
                        if g[2] > max_end:
                            max_end = g[2]
                            longest_gene = g[0]

                    if longest_gene is not None and sj_end > max_end:
                        dog_entry = pd.concat([filtered_genes.loc[longest_gene], filtered_junctions.loc[sj]])
                        #print(f"{sj_start}-{sj_end} count = {filtered_junctions.loc[sj, 4]}")
                        dog_entry.index = ["gene_chrom", "gene_start", "gene_end", "gene_name", "gene_dot", "gene_strand",
                                        "sj_chrom", "sj_start", "sj_end", "sj_name", "sj_count", "sj_strand"]
                        dog_results.append(dog_entry)

                else:  # strand_part == "-"
                    index_start = find_position(gene_coords, sj_end) + 5
                    overlap_genes = []
                    for gene in filtered_genes.index[index_start:0:-1]:
                        gene_start = filtered_genes.loc[gene, "left"]
                        gene_end = filtered_genes.loc[gene, "right"]
                        if gene_start <= sj_end <= gene_end:
                            overlap_genes.append([gene, gene_start, gene_end])
                        elif gene_end < sj_start:
                            break

                    # Find most upstream gene (smallest start)
                    longest_gene = None
                    min_start = float("inf")
                    for g in overlap_genes:
                        if g[1] < min_start:
                            min_start = g[1]
                            longest_gene = g[0]

                    if longest_gene is not None and sj_start < min_start:
                        dog_entry = pd.concat([filtered_genes.loc[longest_gene], filtered_junctions.loc[sj]])
                        #print(f"{sj_start}-{sj_end} count = {filtered_junctions.loc[sj, 4]}")
                        dog_entry.index = ["gene_chrom", "gene_start", "gene_end", "gene_name", "gene_dot", "gene_strand",
                                        "sj_chrom", "sj_start", "sj_end", "sj_name", "sj_count", "sj_strand"]
                        dog_results.append(dog_entry)

        if dog_results:
            all_splice_dogs = pd.DataFrame(dog_results)
            summary_df = summarize_splice_dogs(all_splice_dogs.to_dict(orient="records"), file_label)
            all_summaries.append(summary_df)

            # Export BED file
            sj_only = all_splice_dogs.iloc[:, -6:]
            out_stem = stem.replace("_strand", "_strand_DoG")
            output_path = file_path.with_name(out_stem + ".bed")
            with open(output_path, "w") as f:
                f.write(f"#track name=junctions {strand_part}_strand\n")
                sj_only.to_csv(f, sep="\t", header=False, index=False)

            # Export summary TSV
            summary_df.to_csv(output_path.with_suffix(".tsv"), sep="\t", index=False)
            print(f"Done with {file_label}")
    return all_summaries

#export results
def get_output_path(output_dir, filename):
    return os.path.join(output_dir, filename)

#Merge DoGs that splice out of same gene together and then also have the option to not merge
def merge_splice_dogs_rows(df, file_prefix, output_dir):
    unmerged_cols = [col for col in df.columns if col not in ['gene info', 'Dog_sj_coordinate']]

    #Merge rows that have the same 'gene info' and add up the counts 
    sample_cols = unmerged_cols
    summed = df.groupby("gene info")[sample_cols].sum().reset_index()
    coords = df.groupby("gene info")["Dog_sj_coordinate"].apply(lambda x: ":".join([coord.replace(':+', '').replace(':-', '') for coord in x])).reset_index()

    merged = pd.merge(summed, coords, on="gene info")
    merged = merged[["gene info", "Dog_sj_coordinate"] + sample_cols]

    merged_filename = get_output_path(output_dir, f'merged_{file_prefix}.tsv')
    merged.to_csv(merged_filename, sep='\t', index=False)

    return merged

#normalize splice dogs to TPMs from DESeq and this strand specific 
def normalize_splice_dogs(counts_merged_path, counts_unmerged_path, tpm_path, sample_matrix_path, output_dir, deseq2_dir):
    merged_counts = pd.read_csv(counts_merged_path, sep="\t")
    unmerged_counts = pd.read_csv(counts_unmerged_path, sep="\t")

    #Add +1 to counts
    all_samples = [c for c in merged_counts.columns if c not in ["gene info", "Dog_sj_coordinate"]]
    merged_counts[all_samples] = (merged_counts[all_samples] + 1).astype(int)
    unmerged_counts[all_samples] = (unmerged_counts[all_samples] + 1).astype(int)

    #Load TPMs (keep only TPM columns + gene)
    tpm_df = pd.read_excel(tpm_path)
    tpm_df.rename(columns={"gene": "name"}, inplace=True)
    tpm_cols = [c for c in tpm_df.columns if "(TPMs)" in c]
    tpm_df = tpm_df[["name"] + tpm_cols]

    #Load sample matrix
    sample_matrix = pd.read_csv(sample_matrix_path, sep="\t")
    controls = sample_matrix[sample_matrix["condition"] == "control"]["sample"].tolist()
    tests = sample_matrix[sample_matrix["condition"] != "control"]["sample"].tolist()
    all_samples = controls + tests

    control_tpm_cols = [f"TPM_{s}" for s in controls]
    test_tpm_cols = [f"TPM_{s}" for s in tests]

    #Rename columns for consistency
    tpm_df = tpm_df.rename(columns={f"{s} (TPMs)": f"TPM_{s}" for s in all_samples})

    def process(df, label):
        df["name"] = df["gene info"].str.split(":").str[-1]
        merged = df.merge(tpm_df, on="name", how="left")
        
        expressed = merged.copy().reset_index(drop=True)

        #Remove rows where control TPMs are 0 or NaN (dont need since expressed is filtered)
        expressed = expressed.dropna(subset=control_tpm_cols)
        expressed = expressed[~(expressed[control_tpm_cols] == 0).any(axis=1)]

        #Compute normalized counts
        expressed["mean_TPM"] = expressed[[f"TPM_{s}" for s in all_samples]].mean(axis=1)
        for s in all_samples:
            count_col = s
            tpm_col = f"TPM_{s}"
            norm_col = f"Norm_{s}"
            if count_col in expressed and tpm_col in expressed:
                expressed[norm_col] = (expressed[count_col] / expressed[tpm_col]) * expressed["mean_TPM"]

        #Clean normalized columns
        norm_cols = [f"Norm_{s}" for s in all_samples]
        expressed = expressed[np.isfinite(expressed[norm_cols]).all(axis=1)] #Remove rows where norm_cols = inf from normalizing
        expressed[norm_cols] = expressed[norm_cols].round(0).astype(int) #Round to nearest whole value 
        if "name" in expressed.columns:
            expressed = expressed.drop(columns=['name'])
        
        #Export normalized table
        out_tsv = Path(output_dir) / f"norm.{label}.splice_dogs.tsv"
        expressed.to_csv(out_tsv, sep="\t", index=False)
        print(f"Saved normalized {label} splice DoGs to {out_tsv}")

        #Export DESeq2 input
        deseq_df = expressed[["gene info", "Dog_sj_coordinate"] + norm_cols].copy()
        deseq_df["row_id"] = deseq_df["gene info"] + "|" + deseq_df["Dog_sj_coordinate"]
        deseq_df.set_index("row_id", inplace=True)
        deseq_input_path = Path(output_dir) / f"DESeq2_input_{label}_splice_dogs.tsv"
        deseq_df[norm_cols].to_csv(deseq_input_path, sep="\t")
        print(f"Saved DESeq2 input for {label} splice DoGs to {deseq_input_path}")

        #Export gene info mapping
        mapping_path = Path(output_dir) / f"gene_info_mapping_{label}.tsv"
        deseq_df.reset_index()[["row_id", "gene info", "Dog_sj_coordinate"]].to_csv(mapping_path, sep="\t", index=False)
        print(f"Saved gene info mapping for {label} splice DoGs to {mapping_path}")

        
        #Create deseq2 matrix 
        sample_matrix = pd.read_csv(sample_matrix_path, sep="\t")
        deseq2_matrix = sample_matrix[sample_matrix["sample"].isin([s.replace("Norm_", "") for s in norm_cols])].copy()
        deseq2_matrix["sample"] = "Norm_" + deseq2_matrix["sample"]
        deseq2_matrix = deseq2_matrix.set_index("sample").loc[norm_cols]
        deseq2_matrix_path = Path(deseq2_dir) / "deseq2_matrix.tsv"
        deseq2_matrix.to_csv(deseq2_matrix_path, sep="\t")
        print(f"Saved DESeq2 sample matrix to {deseq2_matrix_path}")
        return expressed


    #Process merged and unmerged splice counts
    merged_result = process(merged_counts, "merged")
    unmerged_result = process(unmerged_counts, "unmerged")

    return merged_result, unmerged_result

def merge_deseq2_with_norm(label, deseq2_output_file, filtered_norm_file, output_dir):
    deseq2_output = pd.read_csv(deseq2_output_file, index_col=0)

    # Extract 'Dog_sj_coordinate' from row_id (everything after the "|")
    deseq2_output["Dog_sj_coordinate"] = deseq2_output.index.str.split("|").str[1]

    # Load normalized counts
    norm_df = pd.read_csv(filtered_norm_file, sep="\t")

    # Merge on Dog_sj_coordinate
    merged = deseq2_output.merge(norm_df, on="Dog_sj_coordinate", how="left")

    # Reorder columns
    ordered_cols = [
        "gene info", "Dog_sj_coordinate",
        "MC1", "MC2", "MC3", "MC4", "MZ1", "MZ3", "MZ4", "MZ5",
        "Norm_MC1", "Norm_MC2", "Norm_MC3", "Norm_MC4",
        "Norm_MZ1", "Norm_MZ3", "Norm_MZ4", "Norm_MZ5",
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    final_cols = [c for c in ordered_cols if c in merged.columns]
    merged = merged[final_cols]

    out_path = Path(output_dir) / f"splice.dogs.{label}.all.results.xlsx"
    merged.to_excel(out_path, engine="openpyxl", index=False, na_rep="NaN")
    print(f"Merging DESeq2 results with normalized counts ({label})...")
    print(f"Saved merged results to {out_path}")
    return merged



def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Splice DoGs pipeline")
    parser.add_argument("project_dir", type=Path, help="Base folder containing input files")
    parser.add_argument("--run_splice", action="store_true", help="Run splice DoG detection")
    parser.add_argument("--run_all", action="store_true", help="Run complete splice DoG pipeline")
    args = parser.parse_args()

    homedir = args.project_dir
    bed_dir = homedir / "bed"
    output_dir = homedir / "splice_dogs_results"
    deseq2_dir = output_dir / "DESeq2_output"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(deseq2_dir, exist_ok=True)

    gtf_candidates = list(homedir.glob("*.annotation.gtf"))
    if not gtf_candidates:
        raise FileNotFoundError(f"No .annotation.gtf file found in {homedir}")
    gtf_path = gtf_candidates[0]

    bed_files = list(bed_dir.glob("*.bed"))
    if not bed_files:
        raise FileNotFoundError(f"No BED files found in {bed_dir}")

    # Sample matrix
    sample_matrix_path = homedir / "sample_matrix.tsv"
    sample_matrix_path = generate_and_validate_sample_matrix(bed_files, sample_matrix_path)

    # Prepare genes
    genes_only = pol2.extract_pol2_genes(gtf_path)
    genes_only = pol2.collapse_isoforms(genes_only)
    genes_only = pol2.remove_embedded_genes(genes_only)
    genes_only = genes_only.sort_values(by=["chromosome", "strand", "left", "right"]).reset_index(drop=True)

    #Remove HLA and IG genes
    hla_ig = genes_only[genes_only["name"].str.startswith("HLA-") | genes_only["name"].str.startswith("IG")]
    genes_only = genes_only.merge(hla_ig, how="left", indicator=True).query('_merge == "left_only"').drop(columns=["_merge"])
    print(f"Loaded {len(genes_only)} cleaned Pol2 genes from {gtf_path}")

    #2 Run splice DoG detection
    if args.run_splice or args.run_all:
        summaries = build_splice_dogs(bed_files, genes_only, output_dir)
        final_summary = pd.concat(summaries, ignore_index=True).drop_duplicates()
        print(f"Done finding Splice DoGs.")

        grouped = final_summary.groupby(["gene info", "Dog_sj_coordinate", "file_name"], as_index=False)["# of sj_count from bed"].first()
        pivoted_summary = grouped.pivot_table(index=["gene info", "Dog_sj_coordinate"], columns="file_name", values="# of sj_count from bed", fill_value=0).reset_index()
        pivoted_summary.columns = [col.replace(".positive:+", ":+").replace(".negative:-", ":-") if isinstance(col, str) else col for col in pivoted_summary.columns]

        sample_matrix = pd.read_csv(sample_matrix_path, sep="\t")
        samples = sample_matrix["sample"].tolist()

        for sample in samples:
            minus_col = f"{sample}:-"
            plus_col = f"{sample}:+"
            unified_col = sample
            pivoted_summary[unified_col] = pivoted_summary.get(minus_col, 0).fillna(0) + pivoted_summary.get(plus_col, 0).fillna(0)

        strand_specific_cols = [f"{s}:{strand}" for s in samples for strand in ["-", "+"]]
        pivoted_summary.drop(columns=[c for c in strand_specific_cols if c in pivoted_summary.columns], inplace=True)

        # Save unmerged splice counts
        unmerged_counts_path = output_dir / "all.splice.dogs.counts.tsv"
        pivoted_summary.to_csv(unmerged_counts_path, sep="\t", index=False)

        # Also save merged splice counts
        merged_df = merge_splice_dogs_rows(df=pivoted_summary, file_prefix="splice.dogs.counts", output_dir=output_dir)
        merged_counts_path = output_dir / "merged_splice.dogs.counts.tsv"
        merged_df.to_csv(merged_counts_path, sep="\t", index=False)

    #3Normalize + DESeq2
    if args.run_all:
        print("\nNormalizing splice DoGs")
        merged_counts_path = output_dir / "merged_splice.dogs.counts.tsv"
        unmerged_counts_path = output_dir / "all.splice.dogs.counts.tsv"
        tpm_path = next(homedir.glob("*.deseq2.results.xlsx"))  # auto-detect TPM file

        merged_result, unmerged_result = normalize_splice_dogs(counts_merged_path=merged_counts_path,
            counts_unmerged_path=unmerged_counts_path, tpm_path=tpm_path, sample_matrix_path=sample_matrix_path,
            output_dir=output_dir, deseq2_dir=deseq2_dir)

        print("\nSplice DoGs normalization complete.")
        print(f"Merged results: {merged_result.shape}, Unmerged results: {unmerged_result.shape}")

        # Paths for DESeq2 inputs
        unmerged_file = output_dir / "DESeq2_input_unmerged_splice_dogs.tsv"
        merged_file = output_dir / "DESeq2_input_merged_splice_dogs.tsv"
        matrix_file = deseq2_dir / "deseq2_matrix.tsv"

        # R script
        r_script = f"""
        setwd("{deseq2_dir}")
        library(DESeq2)

        counts_unmerged <- read.table("{unmerged_file}", header=TRUE, row.names=1, sep="\\t")
        coldata <- read.table("{matrix_file}", header=TRUE, row.names=1, sep="\\t")
        dds_unmerged <- DESeqDataSetFromMatrix(countData = counts_unmerged, colData = coldata, design = ~ condition)
        dds_unmerged <- DESeq(dds_unmerged)
        res_unmerged <- results(dds_unmerged, contrast=c("condition", "test", "control"))
        write.csv(as.data.frame(res_unmerged), file="splice_dogs_deseq2_results.tsv")

        counts_merged <- read.table("{merged_file}", header=TRUE, row.names=1, sep="\\t")
        dds_merged <- DESeqDataSetFromMatrix(countData = counts_merged, colData = coldata, design = ~ condition)
        dds_merged <- DESeq(dds_merged)
        res_merged <- results(dds_merged, contrast=c("condition", "test", "control"))
        write.csv(as.data.frame(res_merged), file="merged_splice_dogs_deseq2_results.tsv")
        """

        r_script_path = deseq2_dir / "run_deseq2.R"
        with open(r_script_path, "w") as f:
            f.write(r_script)
        print(f"R script written to: {r_script_path}")

        # Run DESeq2
        print("Running DESeq2 in R... (this may take a while)")
        subprocess.run(["Rscript", str(r_script_path)], check=True)
        print("DESeq2 finished, results written to DESeq2_output/")

        # Merge R results with file with normalized counts
        merge_deseq2_with_norm(
            label="unmerged",
            deseq2_output_file=output_dir / "DESeq2_output/splice_dogs_deseq2_results.tsv",
            filtered_norm_file=output_dir / "norm.unmerged.splice_dogs.tsv", output_dir=output_dir)

        merge_deseq2_with_norm(
            label="merged",
            deseq2_output_file=output_dir / "DESeq2_output/merged_splice_dogs_deseq2_results.tsv",
            filtered_norm_file=output_dir / "norm.merged.splice_dogs.tsv", output_dir=output_dir)

    end_time = time.time()
    elapsed = end_time - start_time
    minutes, seconds = divmod(int(elapsed), 60)
    print(f"\n Pipeline completed in {minutes} min {seconds} sec")

if __name__ == "__main__":
    main()