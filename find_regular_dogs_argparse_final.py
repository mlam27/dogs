#how to run:
#user must create folder that contains all needed files which are:
# 1) .bw files where folder is named 'bw'
# 2)annotation file where filename ends with 'annotation.gtf'
# 3)tpms from deseq2 where filename ends with '.tpm.xlsx',
# 4)unwanted genes bed file where filename is in format of '*.unwanted.bed'
# 5)dogtools --extract_Pol2_genes.py and gtf.py <-- put these 2 scripts in dogtools folder
# 6)main script which is this script
#User can run script by 'python this.script.py path/to/folder 'coverage' or 'norm_deseq' or 'all'
#user can also use --resume_from after editing sample matrix: example - python this.script.py path/to/folder --resume_from coverage/norm_deseq/all

import argparse
import pandas as pd 
import numpy as np
from pathlib import Path
import pyranges as pr
import glob
import sys
import subprocess
import os
project_dir = Path(__file__).resolve().parent
scripts_dir = project_dir / "dogtools"
from dogtools import extract_Pol2_genes as pol2
from dogtools.gtf import gtf_to_bed

bed_columns = ["chromosome", "left", "right", "name", ".", "strand"]

def define_dog_coordinates(df, plus_end = 100, plus_dog = 1000):
    df = df.copy()
    df['X'] = 0
    df['Y'] = 0
    df.loc[df['strand'] == '+', 'X'] = df['right'] + plus_end
    df.loc[df['strand'] == '+', 'Y'] = df['X'] + plus_dog
    df.loc[df['strand'] == '-', 'Y'] = df['left'] - plus_end
    df.loc[df['strand'] == '-', 'X'] = df['Y'] - plus_dog
    print("DoG coordinates generated.")
    return df

def trim_dogs_by_exons(dog_coordinates, exons_only, dogs_bed, annotation_output_dir):
    dogs_gr = pr.PyRanges(dog_coordinates.rename(columns={
        'chromosome': 'Chromosome', 'left': 'Start', 'right': 'End',
        'name': 'Name', 'strand': 'Strand'}))
    exons_gr = pr.PyRanges(exons_only.rename(columns={
        'chromosome': 'Chromosome', 'left': 'Start', 'right': 'End',
        'name': 'Name', 'strand': 'Strand'}))
    pyranges_output = dogs_gr.join(exons_gr, strandedness="same", suffix="_exon")
    bed_path = annotation_output_dir / "dog_exon_overlap_pyranges.bed"
    pyranges_output.to_bed(bed_path)
    dog_exon_overlap = pd.read_csv(bed_path, sep='\t', header=None, names=[
        'DoG_chr', 'DoG_left', 'DoG_right', 'DoG_name', 'DoG_.', 'DoG_strand',
        'exon_chr', 'exon_left', 'exon_right', 'exon_name', 'exon_.', 'exon_strand'])

     #Plus strand
    plus_strand = dog_exon_overlap[(dog_exon_overlap['DoG_strand'] == '+') & (dog_exon_overlap['exon_strand'] == '+')]
    idx_to_keep = plus_strand.groupby(['DoG_name', 'exon_name'])['exon_left'].idxmin()
    filtered_df = dog_exon_overlap.loc[idx_to_keep].reset_index(drop=True)
    filtered_df = filtered_df.sort_values(by=['DoG_chr', 'DoG_left'], ascending=True).reset_index(drop=True)
    filtered_df_sorted = filtered_df.sort_values(by=['exon_chr','DoG_name', 'exon_left'])
    plus_strand_df = filtered_df_sorted.drop_duplicates(subset='DoG_name', keep='first')
    plus_strand_df = plus_strand_df.sort_values(by=['DoG_chr', 'DoG_left'], ascending=True).reset_index(drop=True)
    trimmed_plus = plus_strand_df.copy()
    trimmed_plus['DoG_right'] = trimmed_plus['exon_left'] - 1

    #Minus strand 
    minus_strand = dog_exon_overlap[(dog_exon_overlap['DoG_strand'] == '-') & (dog_exon_overlap['exon_strand'] == '-')]
    idx_to_keep2 = minus_strand.groupby(['DoG_name', 'exon_name'])['exon_right'].idxmax()
    keep_max = dog_exon_overlap.loc[idx_to_keep2].reset_index(drop=True)
    keep_max = keep_max.sort_values(by=['DoG_chr', 'DoG_left'], ascending=True).reset_index(drop=True)
    keep_max_sorted = keep_max.sort_values(by=['DoG_name', 'exon_right'], ascending=[True, False])
    minus_strand_df = keep_max_sorted.drop_duplicates(subset='DoG_name', keep='first')
    minus_strand_df = minus_strand_df.sort_values(by=['DoG_chr', 'DoG_left'], ascending=True).reset_index(drop=True)
    trimmed_minus = minus_strand_df.copy()
    trimmed_minus['DoG_left'] = trimmed_minus['exon_right'] + 1

    keep = ['DoG_chr', 'DoG_left', 'DoG_right', 'DoG_name', 'DoG_.', 'DoG_strand']
    minus_dog_trimmed = trimmed_minus.copy()
    minus_dog_trimmed = minus_dog_trimmed[keep]
    plus_dog_trimmed = trimmed_plus.copy()
    plus_dog_trimmed = plus_dog_trimmed[keep]

    dog_coords = dogs_bed.copy()
    dog_coords.set_index(['name', 'strand'], inplace=True)
    plus_dog_trimmed.set_index(['DoG_name', 'DoG_strand'], inplace=True)
    minus_dog_trimmed.set_index(['DoG_name', 'DoG_strand'], inplace=True)

    dog_coords.update(plus_dog_trimmed[['DoG_right']].rename(columns={'DoG_right': 'Y'}))
    dog_coords.update(minus_dog_trimmed[['DoG_left']].rename(columns={'DoG_left': 'X'}))
    dog_coords.reset_index(inplace=True)
    plus_dog_trimmed.reset_index(inplace=True)
    plus_dog_trimmed = plus_dog_trimmed[['DoG_chr', 'DoG_left', 'DoG_right', 'DoG_name', 'DoG_.', 'DoG_strand']]
    minus_dog_trimmed.reset_index(inplace=True)
    minus_dog_trimmed = minus_dog_trimmed[['DoG_chr', 'DoG_left', 'DoG_right', 'DoG_name', 'DoG_.', 'DoG_strand']]
    dog_coords[['X', 'Y']] = dog_coords[['X', 'Y']].astype(int)
    dog_coords = dog_coords[['chromosome', 'X', 'Y', 'name', '.', 'strand']] #rearrange

    bad_rows = dog_coords[(dog_coords['X'] > dog_coords['Y']) | (dog_coords['X'] < 0) | (dog_coords['Y'] < 0) | (dog_coords['X'] == dog_coords['Y'])]
    print(f"Number of rows with invalid coordinates: {len(bad_rows)}")
    # Remove those invalid rows
    dog_coords_filtered = dog_coords.drop(bad_rows.index)

    out_path = annotation_output_dir / 'dog_coords_trimmed.bed'
    dog_coords_filtered.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"Done processing file: dog_coords_trimmed.bed")

    return dog_coords_filtered

def remove_unwanted_genes(dog_bed, unwanted_bed_path, output_path):
    if not Path(unwanted_bed_path).exists():
        raise FileNotFoundError(f"ERROR: BED file not found: {unwanted_bed_path}")
    unwanted_df = pd.read_csv(unwanted_bed_path, sep='\t', header=None, names=bed_columns)
    dog_pr = pr.PyRanges(dog_bed.rename(columns={
        'chromosome': 'Chromosome', 'left': 'Start', 'right': 'End',
        'name': 'Name', 'strand': 'Strand'}))
    unwanted_pr = pr.PyRanges(unwanted_df.rename(columns={
        'chromosome': 'Chromosome', 'left': 'Start', 'right': 'End',
        'name': 'Name', 'strand': 'Strand'}))
    overlap = dog_pr.join(unwanted_pr, strandedness="same", suffix="_denovo").df
    to_remove = overlap['Name'].unique()
    final_dogs = dog_bed[~dog_bed['name'].isin(to_remove)]
    final_dogs = final_dogs[final_dogs['left'] != final_dogs['right']]
    final_dogs.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"unwanted genes removed. File exported to: {output_path}")
    return final_dogs

def generate_and_validate_sample_matrix(bw_dir, sample_matrix_path):
    bw_files = glob.glob(os.path.join(bw_dir, "*.bw"))
    print(f"Found {len(bw_files)} negative-strand .bw files.")
    print(f"Found {len(bw_files)} positive-strand .bw files.")
    prefixes = sorted([os.path.basename(f).split('.')[0] for f in bw_files])

    if not os.path.exists(sample_matrix_path):
        print("Sample matrix not found. Generating template...")
        with open(sample_matrix_path, 'w') as f:
            f.write("sample\tcondition\n")
            for p in prefixes:
                f.write(f"{p}\tcontrol_or_test\n")
        raise FileNotFoundError(
            f"Sample matrix made {sample_matrix_path}. Please fill in the 'condition' column (e.g., 'control' or 'test') and rerun the script.")

    df = pd.read_csv(sample_matrix_path, sep='\t')
    missing = [p for p in prefixes if p not in df['sample'].values]
    if missing:
        raise ValueError(f"Sample matrix is missing entries for: {missing}")
    return sample_matrix_path

#Calculate coverage for each sample
def calculate_bw_coverage(bed_path, bw_dir, out_path, sample_matrix_path):
    import pyBigWig
    import glob
    import numpy as np
    import pandas as pd
    import os

    if not os.path.isfile(bed_path):
        raise FileNotFoundError(f"BED file not found at: {bed_path}")
    print(f"Using DoG annotation file: {bed_path}")

    df_matrix = pd.read_csv(sample_matrix_path, sep='\t')
    samples = df_matrix['sample'].tolist()

    neg_bw_files = {os.path.basename(f).split('.')[0]: f for f in glob.glob(os.path.join(bw_dir, "*str1*.bw"))}
    pos_bw_files = {os.path.basename(f).split('.')[0]: f for f in glob.glob(os.path.join(bw_dir, "*str2*.bw"))}

    regions = []
    region_names = []
    with open(bed_path) as bed_handle:
        for line in bed_handle:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                regions.append(fields)
                region_names.append(fields[3])
            else:
                print(f"Skipping line: {line.strip()}")

    results = {name: {} for name in region_names}

    for prefix in samples:
        str1_path = neg_bw_files.get(prefix)
        str2_path = pos_bw_files.get(prefix)
        if not str1_path:
            print(f"Warning: No .str1.bw file found for sample {prefix}")
            continue

        str1 = pyBigWig.open(str1_path)
        str2 = pyBigWig.open(str2_path) if str2_path else None
        stranded = str2 is not None

        for region in regions:
            chrom = region[0]
            start = int(region[1])
            end = int(region[2])
            name = region[3]
            strand = region[5]
            if stranded and strand == '-':
                values = str1.values(chrom, start, end)
            elif stranded and strand == '+':
                values = str2.values(chrom, start, end)
            else:
                values = str1.values(chrom, start, end)

            values = np.array(values, dtype=float)
            if values.size == 0:
                mean = 0
            else:
                values[np.isnan(values)] = 0
                mean = abs(values.mean())

            results[name][prefix] = mean

    df = pd.DataFrame.from_dict(results, orient='index')
    df.index.name = "name"
    df = df.fillna(0).reset_index()
    columns_sorted = ['name'] + sorted([col for col in df.columns if col != 'name'])
    df = df[columns_sorted]

    df.to_csv(out_path, sep="\t", index=False)
    print(f"Results saved to: {out_path}")

#Normalize coverage based off tpms 
def normalize_coverage(coverage_path, tpm_path, sample_matrix_path, out_tsv, deseq_input_path, matrix_out):
    dog_coverage = pd.read_csv(coverage_path, sep='\t')
    deseq_output = pd.read_excel(tpm_path)
    deseq_output.rename(columns={'gene': 'name'}, inplace=True)
    sample_matrix = pd.read_csv(sample_matrix_path, sep='\t')
    controls = sample_matrix[sample_matrix['condition'] == 'control']['sample'].tolist()
    tests = sample_matrix[sample_matrix['condition'] != 'control']['sample'].tolist()

    all_samples = controls + tests
    dog_coverage = dog_coverage.rename(columns={s: f"DoG_{s}" for s in controls + tests})
    deseq_output = deseq_output.rename(columns={f"{s} (TPMs)": f"TPM_{s}" for s in controls + tests})
    merged = dog_coverage.merge(deseq_output, on='name', how='inner')

    dog_cols = [f"DoG_{s}" for s in all_samples]
    tpm_cols = [f"TPM_{s}" for s in all_samples]
    keep_columns = ['name'] + dog_cols + tpm_cols
    merged = merged[keep_columns].copy()

    #Keep only expressed genes (>3)
    expressed = merged[(merged[[f"TPM_{s}" for s in controls]].gt(3).all(axis=1)) & (merged[[f"TPM_{s}" for s in tests]].gt(1).all(axis=1))].copy().reset_index(drop=True)
    expressed['mean_TPM'] = expressed[[f"TPM_{s}" for s in all_samples]].mean(axis=1)

    for s in all_samples:
        dog_col = f"DoG_{s}"
        tpm_col = f"TPM_{s}"
        norm_col = f"{s}_DoG_norm"
        if dog_col in expressed and tpm_col in expressed:
            expressed[norm_col] = (expressed[dog_col] / expressed[tpm_col]) * expressed['mean_TPM']
        else:
            print(f"Skipping {s}: Missing column")

    expressed.to_csv(out_tsv, sep='\t', index=False)
    print(f"Normalized output saved: {out_tsv}")

    #Prepare DESeq2 input
    norm_cols = [f"{s}_DoG_norm" for s in all_samples]
    df = expressed[['name'] + norm_cols].copy()
    df.loc[:, df.columns != 'name'] += 0.1
    df[norm_cols] = (df[norm_cols] * 100).round(0).fillna(0).astype(int)
    df.to_csv(deseq_input_path, sep='\t', index=False)
    print(f"DESeq2 input file saved: {deseq_input_path}")

    matrix = pd.DataFrame({'sample': norm_cols, 'condition': ['control'] * len(controls) + ['test'] * len(tests)})
    matrix.set_index('sample', inplace=True)
    matrix.to_csv(matrix_out, sep='\t')
    print(f"DESeq2 sample matrix saved: {matrix_out}")

    #Final merged output: raw DoG + normalized counts + TPMs
    deseq_input_df = pd.read_csv(deseq_input_path, sep='\t')
    dog_coverage_df = pd.read_csv(coverage_path, sep='\t')
    tpm_clean = deseq_output[['name'] + [f"TPM_{s}" for s in all_samples]].copy()
    merged_df = dog_coverage_df.merge(deseq_input_df, on='name', how='inner')
    merged_df = merged_df.merge(tpm_clean, on='name', how='inner')
    merged_output_path = Path(out_tsv).parent / "norm.rounded.regular.dogs.tsv"
    merged_df.to_csv(merged_output_path, sep='\t', index=False)
    print(f"Final merged output saved to: {merged_output_path}")

def main():
    parser = argparse.ArgumentParser(description="Pipeline for generating and processing regular DoG.")
    parser = argparse.ArgumentParser(description="Pipeline for generating and processing DoG annotations.")
    parser.add_argument("project_dir", type=Path, help="Parent folder containing all required input files and subfolders")
    parser.add_argument("mode", choices=["coverage", "norm_deseq", "all"], help="Which step(s) to run?: 'coverage', 'norm_deseq', or 'all'")
    args = parser.parse_args()

    #Paths to files
    homedir = args.project_dir
    workingdir = homedir / "working_directory"
    bw_dir = homedir / "bw"  # <-- user must place .bw files here

    gtf_candidates = list(homedir.glob("*.annotation.gtf"))
    if not gtf_candidates:
        raise FileNotFoundError(f"No .annotation.gtf file found in {homedir}")
    gtf_path = gtf_candidates[0]

    tpm_candidates = list(homedir.glob("*results.xlsx"))
    if not tpm_candidates:
        raise FileNotFoundError(f"No TPM Excel file (.xlsx) found in {homedir}")
    tpm_path = tpm_candidates[0]

    unwanted = list(homedir.glob("*.unwanted.bed"))
    if not unwanted:
        raise FileNotFoundError(f"No unwanted BED file found in {homedir}")
    unwanted_bed = unwanted[0]

    #Define directories
    output_dir = args.project_dir / "regular_dogs_results"
    annotation_output_dir = output_dir / "annotations"
    deseq2_dir = output_dir / "DESeq2_output"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(annotation_output_dir, exist_ok=True)
    os.makedirs(deseq2_dir, exist_ok=True)

    for f in [gtf_path, tpm_path, unwanted_bed]:
        if not f.exists():
            raise FileNotFoundError(f"Missing required file: {f}")
    if not bw_dir.exists():
        raise FileNotFoundError(f"Missing bigWig directory: {bw_dir}")

    #SAMPLE MATRIX FIRST
    sample_matrix_path = output_dir / 'sample_matrix.tsv'
    sample_matrix_path = generate_and_validate_sample_matrix(bw_dir, sample_matrix_path)

    if args.resume_from in ["dogs", "coverage", "normalization", "deseq2"]:
        #Pol2 gene extraction
        genes_only = pol2.extract_pol2_genes(gtf_path)
        genes_only = pol2.collapse_isoforms(genes_only)
        genes_only = pol2.remove_embedded_genes(genes_only)
        genes_only = genes_only.sort_values(by=["chromosome", "strand", "left", "right"]).reset_index(drop=True)
        out_path = output_dir / "gencode.pol2.filtered.bed"
        genes_only.to_csv(out_path, sep="\t", index=False, header=False)

        print("Gene annotation loaded and filtered for pol2 genes.")

        #Process annotation file, generate dog coordinates, and edit if needed
        dogs = define_dog_coordinates(genes_only)
        dogs_bed = dogs[['chromosome', 'X', 'Y', 'name', '.', 'strand']].copy()
        dogs_bed = dogs_bed[dogs_bed.chromosome != 'chrM'].reset_index(drop=True)
        dogs_bed = dogs_bed.sort_values(by=['chromosome', 'X'], ascending=True).reset_index(drop=True)
        dog_coordinates = dogs_bed.rename(columns={'X': 'left', 'Y': 'right'}).copy()

        exons = gtf_to_bed(gtf_path, feature='exon')
        exons = exons[exons.chromosome != 'chrM'].reset_index(drop=True)
        exons.to_csv(annotation_output_dir / 'exons_only.bed', sep='\t', index=False, header=False)

        find_overlap_and_trim = trim_dogs_by_exons(dog_coordinates, exons, dogs_bed, annotation_output_dir)
        dogs_trimmed = find_overlap_and_trim.rename(columns={'X': 'left', 'Y': 'right'}).copy()
        dogs_trimmed.to_csv(annotation_output_dir / 'dog_coords_trimmed.bed', sep='\t', index=False, header=False)

        final_output = annotation_output_dir / 'dog_coords_without_denovo.bed'
        dog_trimmed_path = annotation_output_dir / 'dog_coords_trimmed.bed'
        dog_df = pd.read_csv(dog_trimmed_path, sep='\t', header=None, names=bed_columns)
        dog_no_denovo = remove_unwanted_genes(dog_df, unwanted_bed, final_output)

        print("Done prepping dog coordinates!")
    else:
        final_output = annotation_output_dir / 'dog_coords_without_denovo.bed'

    #Run coverage function
    out_path = output_dir / 'dog_coverage.tsv'
    if args.mode in ['coverage', 'all']:
        calculate_bw_coverage(final_output, bw_dir, out_path, sample_matrix_path)

    #Run normalization and DESeq2
    if args.mode in ["norm_deseq", "all"]:
        normalize_coverage(
            coverage_path=out_path,
            tpm_path=tpm_path,
            sample_matrix_path=sample_matrix_path,
            out_tsv=deseq2_dir / 'norm.raw.regular.dogs.tsv',
            deseq_input_path=deseq2_dir / 'DESeq2_input.tsv',
            matrix_out=deseq2_dir / 'dog_matrix.tsv')

        if not r_script_path.exists():
            r_script_contents = """\
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
matrix_file <- args[2]
output_csv <- args[3]

library(DESeq2)

counts <- read.table(input_file, header=TRUE, row.names=1, sep="\\t")
coldata <- read.table(matrix_file, header=TRUE, row.names=1, sep="\\t")

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "test", "control"))

write.csv(as.data.frame(res), file=output_csv)
"""
        r_script_path.write_text(r_script_contents.strip())
        print(f"Saved DESeq2 R script to: {r_script_path}")

    #Run R script
        try:
            subprocess.run([
                "Rscript",
                str(r_script_path),
                str(counts_path),
                str(matrix_path),
                str(output_csv)
            ], check=True)
            print(f"DESeq2 analysis completed. Results saved to: {output_csv}")
        except subprocess.CalledProcessError as e:
            print(f"Rscript failed with exit code {e.returncode}")
            raise

        #Merge results back in Python
        file_to_merge = deseq2_dir / 'norm.rounded.regular.dogs.tsv'
        if output_csv.exists():
            try:
                dog = pd.read_csv(file_to_merge, sep='\t')
                deseq2_output = pd.read_csv(output_csv)
                if deseq2_output.columns[0] != "name":
                    deseq2_output = deseq2_output.rename(columns={deseq2_output.columns[0]: "name"})
                merged = dog.merge(deseq2_output, on='name', how='left')

                merged.to_excel(final_output_excel, engine="openpyxl", index=False)
                print(f"Merged DESeq2 results saved to: {final_output_excel}")
            except Exception as e:
                print(f"Failed to merge or write final Excel file: {e}")
        else:
            print("DESeq2 output CSV not found. Cannot proceed with merging.")

if __name__ == "__main__":
    main()