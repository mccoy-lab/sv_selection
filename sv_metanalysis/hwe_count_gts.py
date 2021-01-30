#!/usr/bin/env python3

"""
Usage: ./hwe_count_genotypes.py <VCF> <manifest file>

Generates files with counts of samples that are hom ref, hom alt, and het for each SV,
within each population. Used as input into R's HardyWeinberg package. Autosomal and
X chromosome variants are treated separately. Y chromosome variants are ignored.
"""

import sys
from cyvcf2 import VCF
import pandas as pd


"""
Generates autosomal genotype counts for a list of samples.
(i.e. samples that belong to one population)
"""
def count_gts(sample_list):
    # load VCF and choose only samples of interest
    vcf = VCF(sys.argv[1])
    vcf.set_samples(sample_list)
    
    # initialize lists
    ids = []
    hets = []
    refs = []
    alts = []
    missing = []

    # count # of samples of each genotype for every variant
    for variant in vcf:
        if (variant.CHROM != "chrX") and (variant.CHROM != "chrY"): # skip X, Y chromosome variants
            ids.append(variant.ID)
            hets.append(variant.num_het)
            refs.append(variant.num_hom_ref)
            alts.append(variant.num_hom_alt)
            missing.append(variant.num_unknown)

    # create dataframes with gt counts
    df = pd.DataFrame(
        {"SV" : ids,
        "ref" : refs,
        "het" : hets,
        "alt" : alts,
        "missing" : missing
        })
    return(df)

"""
Generates X chromosome genotype counts for a list of samples.
(i.e. samples that belong to one population)
Only outputs ref/alt genotypes for male samples, since they're X chromosome haploids.
"""
def count_gts_xchr(sample_list, sex):
    # load VCF and choose only samples of interest
    vcf = VCF(sys.argv[1])
    vcf.set_samples(sample_list)
    
    # initialize lists
    ids = []
    hets = []
    refs = []
    alts = []
    missing = []
        
    # count # of samples of each genotype for every variant
    for variant in vcf:
        if variant.CHROM == "chrX": # only look at variants on X chromosome
            ids.append(variant.ID)
            hets.append(variant.num_het)
            refs.append(variant.num_hom_ref)
            alts.append(variant.num_hom_alt)
            missing.append(variant.num_unknown)

    # create dataframes with gt counts
    df = pd.DataFrame(
        {"SV" : ids,
        "ref" : refs,
        "het" : hets,
        "alt" : alts,
        "missing" : missing
        })
    if sex == "male":
        # incorporate het samples into alt allele calls
        df["alt"] = df["het"] + df["alt"]
        return(df[["SV","ref","alt","missing"]])
    else:
        return(df[["ref","het","alt","missing"]])


"""
Create a dataframe of samples and their corresponding populations.
"""

samples = pd.read_table(sys.argv[2], header=None, names=["ident","cov","sex","pop"], index_col=4)
# split sample path column to return only the sample ID
def split_col(ident):
    return((ident.split(".")[0]).split("/")[1])
samples["ident"] = samples["ident"].apply(lambda x: split_col(x))

populations = samples["pop"].unique()


"""
Loop through populations, counting autosome and X chromosome genotypes for each.
Output genotype counts as separate files.
"""

output_dir = "gt_counts/"

for pop in populations:
    # get IDs for samples in population
    roi = samples.loc[:,"pop"] == pop
    subset = samples.loc[roi,:]
    pop_ids = subset["ident"].tolist()
    
    # get IDs for male and female samples in population
    m_roi = subset.loc[:,"sex"] == "male"
    m_subset = subset.loc[m_roi,:] # male samples
    m_pop_ids = m_subset["ident"].tolist()
    f_subset = subset.loc[-m_roi,:] # female samples
    f_pop_ids = f_subset["ident"].tolist()
    
    # get genotype counts for autosomal variants
    autosome_counts = count_gts(pop_ids)
    autosome_counts.to_csv(output_dir + pop + "_gtcounts_autosome.txt", sep="\t", index=False, header=False)
    
    # get genotype counts for X chromosome variants
    xchr_m_counts = count_gts_xchr(m_pop_ids, "male")
    xchr_f_counts = count_gts_xchr(f_pop_ids, "female")
    xchr_counts = pd.concat([xchr_m_counts, xchr_f_counts], axis=1)
    # add counts of missing GT calls from male and female samples
    xchr_counts.iloc[:,7] = xchr_counts.iloc[:,3] + xchr_counts.iloc[:,7]
    xchr_counts.iloc[:,[0,1,2,4,5,6,7]].to_csv(output_dir + pop + "_gtcounts_xchr.txt", sep="\t", index=False, header=False)
