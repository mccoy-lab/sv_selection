#!/usr/bin/env python3

"""
Usage: ./calculate_sgdp_af_02.py

Calculates ref and alt allele frequencies for an SGDP VCF using cyvcf2.
"""

import sys
import pandas as pd
from cyvcf2 import VCF

# sgdp metadata
# converted to a CSV using the first sheet in: https://static-content.springer.com/esm/art%3A10.1038%2Fnature18964/MediaObjects/41586_2016_BFnature18964_MOESM205_ESM.xlsx
meta = pd.read_csv(open("sgdp_meta.csv"))
meta_subset = meta.iloc[0:300,] # get rid of empty rows

# get names of unique populations in dataset
populations = set(meta_subset.loc[:, "Population ID"].tolist())

# create output file with population-specific AFs
w = open("sgdp_afs_population.txt", "w")
w.write("Variant\tAF\tident\tPopulation\tRegion\n") # headers

# iterate through populations to calculate population-specific AF
for pop in populations:
    roi = (meta_subset["Population ID"] == pop)
    sample_ids = meta_subset.loc[roi, "Sample ID (SGDP)"].tolist()
    region = meta_subset.loc[roi, "Region"].tolist()[0]

    vcf =  VCF("sgdp_chr14_IDs.vcf.gz", samples = sample_ids)
    for variant in vcf:
        ident = variant.ID
        aaf = str(variant.aaf) # alt allele frequency
        raf = str(1 - variant.aaf) # ref allele frequency
        # write to output file
        w.write(ident + "\t" + aaf + "\t" + "ALT" + "\t" + pop + "\t" + region + "\n")
        w.write(ident + "\t" + raf + "\t" + "REF" + "\t" + pop + "\t" + region + "\n")
    vcf.close()