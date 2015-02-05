#!/usr/bin/env python

# About:
# This script prints a one line header close to what the old SAVI produced

# Usage:
# arg1 is the comma-delimited list of sample names (IN THE ORDER IN WHICH YOU RAN THEM IN SAVI)
# arg2 is the comma-delimited list of comparisons (which are themselves colon delimited) (IN THE ORDER IN WHICH YOU RAN THEM IN SAVI) 

# E.g., 
# make_header.py 1,2,3 2:1,3:1

import sys

namelist = sys.argv[1].split(",")
complist = sys.argv[2].split(",")

# print first field
print('#chrom'),

# static fields
mylist1=['pos', 'id', 'ref', 'var', 'AA', 'CDS', 'cosmic_number_of_samples', 'Effect', 'Effect_Impact', 'Functional_Class', 'Codon_Change', 'Amino_Acid_Change', 'Amino_Acid_length', 'Gene_Name', 'Transcript_BioType', 'Gene_Coding', 'Transcript_ID', 'Exon_Rank', 'Genotype_Number', 'SnpEff_errors', 'SnpEff_warnings', 'meganormal_id', 'number_of_mutations_per_meganormal_id']

# print static fields
for i in mylist1:
	print("\t" + i),

# print names in dictionary order
for i in sorted(namelist):
	print("\t" + i + "f"),
	print("\t" + i + "f_lower"),
	print("\t" + i + "f_upper"),

for i in complist:
	elt1, elt2 = i.split(":")
	print("\t" + elt1 + "f-" + elt2 + "f"),
	print("\t" + elt1 + "f-" + elt2 + "f_lower"),
	print("\t" + elt1 + "f-" + elt2 + "f_upper"),

# more static fields
mylist2=['reference_SNP_id(RS)']
for i in mylist2:
	print("\t" + i),

# print names in dictionary order
for i in sorted(namelist):
	print("\t" + i + "_presence_bool"),
	print("\t" + i + "_presence_posterior"),

# more static fields
mylist3=['cbio_somatic_verified_mutation_count', 'meganormal_186_TCGA_src', 'strand', 'cbio_total_mutation_count', 'dbSNPBuild_for_RS']
for i in mylist3:
	print("\t" + i),

# print samples in the order in which they were run in savi
for i in namelist:
	print("\t" + i + "_totdepth"),
	print("\t" + i + "_refdepth"),
	print("\t" + i + "_vardepth"),
	print("\t" + i + "_ref_ave_qual"),
	print("\t" + i + "_var_ave_qual"),
	print("\t" + i + "_ref_forward_depth"),
	print("\t" + i + "_ref_reverse_depth"),
	print("\t" + i + "_var_forward_depth"),
	print("\t" + i + "_var_reverse_depth"),
