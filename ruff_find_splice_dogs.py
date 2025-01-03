import pandas as pd
from pathlib import Path
from gtf import gtf_to_bed
import glob
import os
import argparse
import numpy as np

files = glob.glob('/Users/meggielam/Desktop/ruff_2024/practice_overlap/*')
min_counts = 10
min_bases = 100


### Where you keep this script (and associated files) ###
homedir = Path("/Users/meggielam/Desktop/ruff_2024")

scriptdir = homedir / "ruff_2024"
annotationdir = homedir / "annotations"
gtf_path = annotationdir / 'gencode.v41.annotation.gtf'
#gtf_path = annotationdir / 'hg38gencode.v42_plus_Akata_inverted.gtf'
bed_path = gtf_path.parent / (gtf_path.stem + '.bed')

### The script converts GTF to bed unless that has already been done by previously running this script ###
if bed_path.exists(): #If bed file exists, read in the bed file using pandas library (pd=pandas)and place it in a dataframe called "genes"#
    genes = pd.read_table(bed_path)
else: #if bed file isn't there, convert the gtf to bed and loop back to 'if' statement#
    genes = gtf_to_bed(gtf_path, feature='transcript')

### Sort gene coordinates based off start from low to high ###
genes = genes.sort_values(by=['left'], ascending=True).reset_index(drop=True)

### Take most left and most right value for each gene ###
#If gene has isoforms, take min 'left' value and max 'right' value and make 1 row with that gene and those coordinates only
grouped_genes = genes.groupby('name').agg({
    'chromosome': 'first',   # Keep the first occurrence (assuming it's the same for all rows of the same name)
    'left': 'min',           # Take the minimum 'left' value
    'right': 'max',          # Take the maximum 'right' value
    'type': 'first',         # Keep the first occurrence
    'strand': 'first'        # Keep the first occurrence (assuming it's the same for all rows of the same name)
}).reset_index()

#Resort gene coordinates based off start position from low to high before removing embedded genes
genes = genes.sort_values(by=['left'], ascending=True).reset_index(drop=True)


###Remove embedded genes(Both strands)###
def update_gene_range(current_gene):
    return current_gene['left'], current_gene['right'] 

embedded_genes = []
gene_left = genes.iloc[0]['left']
gene_right = genes.iloc[0]['right']
reference_chromosome = genes.iloc[0]['chromosome']
genes.append(genes.iloc[0])

i = 0
while i < len(genes):
    current_gene = genes.iloc[i]
    if current_gene['strand'] == '+': #For positive strand#
        if (current_gene['left'] >= gene_left and current_gene['right'] <= gene_right and current_gene['chromosome'] == reference_chromosome):
            print("Condition met for '+': Adding to embedded_genes and updating 'genes'")
            embedded_genes.append(current_gene)
            genes = genes.drop(genes.index[i]).reset_index(drop=True)
            continue
    elif current_gene['strand'] == '-': #For negative strand#
        if (current_gene['left'] <= gene_left and current_gene['right'] >= gene_right and current_gene['chromosome'] == reference_chromosome):
            print("Condition met for '-': Adding to embedded_genes and updating 'genes'")
            embedded_genes.append(current_gene)
            genes = genes.drop(genes.index[i]).reset_index(drop=True)
            continue
    
    #Condition not met: Updating gene range 
    gene_left, gene_right = update_gene_range(current_gene)
    reference_chromosome = current_gene['chromosome']
    i += 1

embedded_genes_df = pd.DataFrame(embedded_genes)



### The point of this 'def' function is to make the entire script faster so the loop does not have to go through every line in the annotation file ###
### This function will break the annotation file into sections - lower, middle, and upper - ###
### If the splice junction from the bed file is not in the lower_bound, then it will change the conditions of the lower_bound to continue the loop ###
def find_position(arr, val):
    '''Finds the position of a junction in an array of coordinates. Must be a sorted array'''

    lower_bound = 0
    upper_bound = len(arr)
    while lower_bound < upper_bound:
        middle = (lower_bound + upper_bound) // 2 #defines middle section; '//2' = rounds down to nearest whole number#
        if val >= arr[middle]: #Looks for splice junction in middle section of annotation file, if splice junction is not within that section, update boundaries and continue search#
            lower_bound = middle + 1
        else:
            upper_bound = middle
    return lower_bound - 1 #If the lower_bound is no longer < upper_bound, then go start over with new boundaries#

### This 'def' function will grab splice dogs from positive strand bed file ###
def get_dog_positive(junctions, genes, chromosome):
    '''Output positive strand DoGs'''
    ## Go through bed file and grab splice junctions from positive strands ##
    dog_genes = [] #Empty dataFrame for dog genes to go in
    genes = genes[(genes["chromosome"] == chromosome) & (genes["strand"] == '+')].sort_values('left') #Grabs columns "chromosome", "strand", and sorts the genes based off start position#
    junctions = junctions[(junctions[0] == chromosome) & (junctions[4] >= min_counts)] #column 4 has the counts; Grabs junctions that are >=min count (min count = 0)#
    genes.index = range(len(genes.index))
## Loop through 'genes' to identify splice dogs based on boundries that are set ##
    gene_lefts = list(genes['left']) #Take 'start' position of genes ##
    ## Get start and end positions of each splice junction from bed file and call it sj_start/end ##
    for sj in junctions.index:
        sj_left = junctions.loc[sj, 1]
        sj_end = junctions.loc[sj, 2]

        index_left = find_position(gene_lefts, sj_left)
        for gene in genes.index[index_left:]:
            gene_left = genes.loc[gene, "left"]
            gene_end = genes.loc[gene, "right"]
            if gene_left <= sj_left <= gene_end:
                if sj_end > (gene_end + args.min_bases):
                    dog_genes.append(pd.concat([genes.loc[gene], junctions.loc[sj]]))
            elif gene_left > sj_end:
                break

    return pd.DataFrame(dog_genes)

def get_dog_negative(junctions, genes, chromosome):

    '''Output negative strand DoGs'''
    dog_genes = []
    genes = genes[(genes["chromosome"] == chromosome) & (genes["strand"] == '-')].sort_values('right')
    junctions = junctions[(junctions[0] == chromosome) & (junctions[4] >= args.min_counts)]
    genes.index = range(len(genes.index))

    gene_lefts = list(genes['right'])
    for sj in junctions.index:
        sj_left = junctions.loc[sj, 1]
        sj_end = junctions.loc[sj, 2]

        # Double check this..finds the first position then moves back a few rows and cycles through.
        index_left = find_position(gene_lefts, sj_end) + 5

        # Cycle through rows backwards
        for gene in genes.index[index_left:0:-1]:
            gene_left = genes.loc[gene, "left"]
            gene_end = genes.loc[gene, "right"]
            if gene_left <= sj_end <= gene_end:
                if sj_left < (gene_start - min_bases):
                    dog_genes.append(pd.concat([genes.loc[gene], junctions.loc[sj]]))
            elif gene_end < sj_left:
                break

    return pd.DataFrame(dog_genes)


chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['X','Y']
files = glob.glob(args.files)

# SJ bed file names must have the word "negative" or "positive" in them to distinguish strand. i.e. 'MZ.positive_strand.bed'
for sj_bedfile in files:
    print(f"Processing file: {sj_bedfile}")
    if 'negative' in sj_bedfile:
        negsj_path = Path(sj_bedfile)
        sjs_neg = pd.read_table(negsj_path, comment='#', header=None)
        neg = pd.concat([get_dog_negative(sjs_neg, genes, chromosome) for chromosome in chromosomes])
        neg[3] = neg['name']
        negpath_out = negsj_path.parent / (negsj_path.stem + '.dogs_only.bed')
        with open(negpath_out, 'w') as outfile:
            outfile.write("#track name=junctions negative_strand\n")
            neg.iloc[:, 6:].to_csv(outfile, sep='\t', header=None, index=None)

    elif 'positive' in sj_bedfile:
        possj_path = Path(sj_bedfile)
        sjs = pd.read_table(possj_path, comment='#', header=None)
        pos = pd.concat([get_dog_positive(sjs, genes, chromosome) for chromosome in chromosomes])
        pos[3] = pos['name']
        pospath_out = possj_path.parent / (possj_path.stem + '.dogs_only.bed')
        with open(pospath_out, 'w') as outfile:
            outfile.write("#track name=junctions positive_strand\n")
            column_labels = ['1', '2', '3', '4', '5', '6']
            outfile.write('\t'.join(column_labels) + '\n')
            pos.iloc[:, 6:].to_csv(outfile, sep='\t', header=None, index=None)
    print("~~~~DONE~~~~")









# ##to normalize to gene expression
# import pandas as pd 
# import glob

# gene_expression = "/Users/meggielam/Desktop/ruff_2024/mutu_polyA_stranded_facs_bed/zta.summary.genes.tpm.txt"
# zta_genes = pd.read_table(gene_expression, index_col=0)

# #to go through all files that start with 'MC' 
# directory_path = "/Users/meggielam/Desktop/ruff_2024/mutu_polyA_stranded_facs_bed/"
# mutu_dog_files = glob.glob(directory_path+"M*.*.dogs_only.bed")

# for file in mutu_dog_files: #mutu dog only files 
#     name = file.split('.')[0].split('/')[-1] #will take sample name from filename 
#     df = pd.read_csv(file, sep="\t",comment='#', index_col=3,header=None) #reads the files
#     df = pd.DataFrame(df[4]) # Number of DoG junctions
#     df['expression'] = zta_genes[name] #makes a column called expression in df and then indexes thru zta_genes file by 'name'
#     df['normalized_dog']= df[4] / df['expression'] #normalizes by gene expression number from zta_genes file
#     df.columns = ["dogs", "expression", "normalized_dogs"] #rename columns
#     df.to_csv(file+'.norm.tsv', sep='\t')  #saves and adds 'norm.tsv' to file name


