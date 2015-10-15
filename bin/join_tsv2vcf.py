#!/usr/bin/env python

# About:
# hacky script join dna tsv and rna vcf on common variants

# Usage:
# join_tsv2vcf.py dna.tsv rna.vcf > out.tsv

import sys, re, zipfile

### Variables
d = {}			# dict with keys as variants, values as rna depths 
num_samples = 0 	# num of samples in vcf file
sdpindex = 0		# index of 'SDP'
adindex = 0		# index of 'AD'
isfirst = 1		# boolean denoting first line

### Main

# loop thro rna vcf
with open(sys.argv[2], "r") as f:
# with open(zipfile.ZipFile(sys.argv[2]), "r") as f:
	for line in f:
		# ignore comments
		if ( line[0:1] != "#" ):

			# print(line),

			# first line
			if (isfirst):
				# get list of elts in format field
				num_samples = len(line.split()) - 9
				
				# get indices
				sdpindex = line.split()[8].split(':').index('SDP')
				adindex = line.split()[8].split(':').index('AD')
	 
				# turn off flag
				isfirst = 0

			# create key
			mykey = ":".join(map(line.split().__getitem__, [0,1,3,4]))

			# string for rna depths
			rnadepths = ""
			for i in range(9,len(line.split())):
				# extract total depth, variant depth
				rnadepths += line.split()[i].split(':')[sdpindex] + "\t" + line.split()[i].split(':')[adindex] + "\t"
			d[mykey] = rnadepths.rstrip()
				
# loop thro dna tsv
with open(sys.argv[1], "r") as f:
	for line in f:
		# print double # header
		if ( line[0:2] == "##" ):
			print(line),
		# if single # header, augment with rna header
		elif ( line[0:1] == "#" ):
			print(line.rstrip() + "\t" + "\t".join(["rna_refdepth_" + str(i) + "\trna_altdepth_" + str(i) for i in range(1,num_samples+1)]))
		else:
			mykey2 = ":".join(line.split()[0:4])
			if mykey2 in d:
				print(line.rstrip() + "\t" + d[mykey2])
			else:
				print(line.rstrip() + "\t" + "\t".join(["-\t-" for i in range(1,num_samples+1)]))
