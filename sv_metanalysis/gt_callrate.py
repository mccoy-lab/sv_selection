#!/usr/bin/env python3

"""
Usage: ./gt_callrate.py <VCF> <manifest file>

Determines genotyping rate for SVs.
For variants on the Y chromosome, only male samples are considered.
"""

import sys
from cyvcf2 import VCF
import pandas as pd


# for non-Y chromosome variants
def get_call_rate_nonY(sample_list):
    vcf = VCF(path_to_vcf)
    vcf.set_samples(sample_list)

    # initialize lists
    ids = []
    call_rate = []
    
    # count # of samples of each genotype for every variant
    for variant in vcf:
        if variant.CHROM != "chrY": # ignore variants on Y chromosome
            ids.append(variant.ID)
            call_rate.append(variant.call_rate)
    
    # create dataframe
    df = pd.DataFrame(
        {"SV" : ids,
        "call_rate" : call_rate
        })
    return(df)

# for Y chromosome variants
def get_call_rate_Ychr(sample_list):
    vcf = VCF(path_to_vcf)
    vcf.set_samples(sample_list)

    # initialize lists
    ids = []
    call_rate = []

    # count # of samples of each genotype for every variant
    for variant in vcf:
        if variant.CHROM == "chrY": # only look at variants on Y chromosome
            ids.append(variant.ID)
            call_rate.append(variant.call_rate)
    
    # create dataframe
    df = pd.DataFrame(
        {"SV" : ids,
        "call_rate" : call_rate
        })
    return(df)


"""
Apply functions to VCF.
"""
    
path_to_vcf = sys.argv[1] # VCF file
path_to_manifest = sys.argv[2] # Paragraph manifest file

# create lists of samples from manifest file
samples = pd.read_table(path_to_manifest, header=None, names=["ident","cov","sex","pop"], index_col=4)

# split sample path column to return only the sample ID
def split_col(ident):
    return((ident.split(".")[0]).split("/")[1])
samples["ident"] = samples["ident"].apply(lambda x: split_col(x))

all_pop_ids = samples["ident"].tolist() # all sample IDs
m_roi = samples.loc[:,"sex"] == "male" # male samples only
m_subset = samples.loc[m_roi,:]
m_pop_ids = m_subset["ident"].tolist() 

# get call rates
nonY_callrate = get_call_rate_nonY(all_pop_ids)
Ychr_callrate = get_call_rate_Ychr(m_pop_ids)

# combine all variants into one output file
all_callrates = pd.concat([nonY_callrate, Ychr_callrate])
all_callrates.to_csv("genotyping_callrates.txt", sep="\t", index=False, header=True)