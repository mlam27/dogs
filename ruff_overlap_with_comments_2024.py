import pandas as pd
from pathlib import Path
from gtf import gtf_to_bed
import glob
import os
import argparse
import numpy as np



# Path to both negative and positive splice junction bed files
files = glob.glob('/Volumes/FlemingtonLab/lab_experiments/1_EBV_reactivation_project/1_RNA_seq/50_Encode_RBP_shRNA/rMATS-Cntl_Test_erik_runs/1_HepG2/*/*/1_STAR/2_splicing/*negative*')


# Minimum number of splice junction counts required
min_counts = 0

# Where you keep this script (and associated files)
homedir = Path("/Volumes/Flemington_Lab_Documents/21_lab_member_directories/10_Meggie/")
scriptdir = homedir / "ruff_2024/"
annotationdir = scriptdir / "annotations"
gtf_path = annotationdir / 'gencode.v41.annotation.gtf'
bed_path = gtf_path.parent / (gtf_path.stem + '.bed')

### The script converts GTF to bed unless that has already been done by previously running this script ###
### If bed file exists, read in the bed file using pandas library (pd=pandas)and place it in a dataframe called "genes" ###
if bed_path.exists(): 
    genes = pd.read_table(bed_path)
else: #if bed file isn't there, convert the gtf to bed and loop back to 'if' statement
    genes = gtf_to_bed(gtf_path, feature='gene')

### The point of this 'def' function is to make the entire script faster so the loop does not have to go through every line in the annotation file ###
## This function will break the annotation file into sections - lower, middle, and upper -
## If the splice junction from the bed file is not in the lower_bound, then it will change the conditions of the lower_bound to continue the loop ##
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


### This 'def' function will grab splice dogs from positive strand bed file ###
def get_dog_positive(junctions, genes, chromosome):
    '''Output positive strand DoGs'''
    ## Go through bed file and grab splice junctions from positive strands ##
    dog_genes = [] #empty dataFrame for dog genes to go in
    genes = genes[(genes["chromosome"] == chromosome) & (genes["strand"] == '+')].sort_values('left') #grabs columns "chromosome", "strand", and sorts the genes based off start position
    junctions = junctions[(junctions[0] == chromosome) & (junctions[4] >= min_counts)] #grabs junctions that are >=min count (min count = 0)
    genes.index = range(len(genes.index)) 

    ## Loop through 'genes' to identify splice dogs based on boundries that are set ##
    gene_starts = list(genes['left']) #Take 'start' position of genes
    ## Get start and end positions of each splice junction from bed file and call it sj_start/end ##
    for sj in junctions.index: 
        count_it = False #track overlap event
        sj_start = junctions.loc[sj, 1] #set variable for splice junction start
        sj_end = junctions.loc[sj, 2] #set variable for splice junction end 

        index_start = find_position(gene_starts, sj_start) #Use 'find_position' function to find overlapping genes
        overlap_genes = [] #Puts overlap genes into separate list
        #Find overlap genes, if it doesn't overlap then break loop. If overlap genes are found, put into 'overlap_genes' list
        for gene in genes.index[index_start:]:
            gene_start = genes.loc[gene, "left"]
            gene_end = genes.loc[gene, "right"]
            if gene_start <= sj_start <= gene_end:
                count_it = True
                '''store genes, gene_start, gene_end in a list 
                [ 
                    [TP53, 150002, 200000],
                    [CDKN1A, 160345, 653545],
                ]
                '''
                overlap_genes.append([gene, gene_start, gene_end])
            elif gene_start > sj_end:
                break
        #Go through overlap_genes list and find the overlap_gene that is farthest downstream 
        max_downstream_end = 0
        for overlap_gene in overlap_genes:
            if overlap_gene[2] > max_downstream_end:
                max_downstream_end = overlap_gene[2]
                longest_gene = overlap_gene[0]
        #If the overlap_gene is the farthest downstream, then put it into 'dog_genes' dataFrame
        if count_it is True and sj_end > max_downstream_end:
            dog_genes.append(pd.concat([genes.loc[longest_gene], junctions.loc[sj]]))



    return pd.DataFrame(dog_genes)

### Find splice dogs, but in the negative strand ###
## Did not consider overlap genes in negative strand yet since strategy for positive strand did not work ##
def get_dog_negative(junctions, genes, chromosome):

    '''Output negative strand DoGs'''
    dog_genes = []
    genes = genes[(genes["chromosome"] == chromosome) & (genes["strand"] == '-')].sort_values('right')
    junctions = junctions[(junctions[0] == chromosome) & (junctions[4] >= min_counts)]
    genes.index = range(len(genes.index))

    gene_starts = list(genes['right'])
    for sj in junctions.index:
        sj_start = junctions.loc[sj, 1]
        sj_end = junctions.loc[sj, 2]

        # Double check this..finds the first position then moves back a few rows and cycles through.
        index_start = find_position(gene_starts, sj_end) + 5

        # Cycle through rows backwards
        for gene in genes.index[index_start:0:-1]:
            gene_start = genes.loc[gene, "left"]
            gene_end = genes.loc[gene, "right"]
            if gene_start <= sj_end <= gene_end:
                if sj_start < gene_start:
                    dog_genes.append(pd.concat([genes.loc[gene], junctions.loc[sj]]))
            elif gene_end < sj_start:
                break

    return pd.DataFrame(dog_genes)



#### Goal output file format ####

# chromosome_start_stop_strand_geneName [tab] count
# chr1_10000_10200_-_TP53    102



# This is designed to work only for human samples so far - just forms a list of chromosomes
chromosomes = ['chr' + str(i) for i in range(22)] + ['X','Y']

# SJ bed file names must have the word "negative" or "positive" in them to distinguish strand. i.e. 'MZ.positive_strand.bed'
for sj_bedfile in files:
    rnabp = sj_bedfile.split('/')[-1].split('.')[0].split('_')[1]
    output_filename = rnabp + ".dogs.tsv"
    
    
    

    negsj_path = Path(sj_bedfile)
    sjs_neg = pd.read_table(negsj_path, comment='#', header=None)
    neg = pd.concat([get_dog_negative(sjs_neg, genes, chromosome) for chromosome in chromosomes])
    neg[3] = neg['name']
    neg = pd.DataFrame(neg,dtype=str)
    neg = neg.reset_index()
    l = []
    for i in neg.index:
        string = '_'.join([neg.loc[i,0], neg.loc[i,1], neg.loc[i,2], neg.loc[i,5], neg.loc[i,"name"]])
        l.append(string)
    neg.index = l 
         
    
    possj_path = Path(sj_bedfile.replace('negative', 'positive'))
    sjs = pd.read_table(possj_path, comment='#', header=None)
    pos = pd.concat([get_dog_positive(sjs, genes, chromosome) for chromosome in chromosomes])
    pos[3] = pos['name']
    pos = pd.DataFrame(pos,dtype=str)
    pos = pos.reset_index()
    l = []
    for i in pos.index:
        string = '_'.join([pos.loc[i,0], pos.loc[i,1], pos.loc[i,2], pos.loc[i,5] ,pos.loc[i,"name"]])
        l.append(string)
    pos.index= l 
    all_junctions = list(neg.index) + list(pos.index)
    new = pd.DataFrame(index=all_junctions)
    neg = neg.reset_index().drop_duplicates('level_0').set_index('level_0')
    pos = pos.reset_index().drop_duplicates('level_0').set_index('level_0')
    
    new = pd.concat([neg[4], pos[4]])
    pos.to_excel("new.positive_only."+ output_filename +'.xlsx')
    print(output_filename, "DONE")












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


