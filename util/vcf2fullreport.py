#!/usr/bin/env python

# About:
# This script takes a vcf file as input, piped in, and converts it to a tab-delimited text file, 
# preserving all information

# Usage:
# cat myfile.vcf | vcf2report.py > myfile.txt
# To suppress lines beginning with "##" (i.e., the vcf header), run:
# cat myfile.vcf | vcf2report.py 0 > myfile.txt

import sys
import re

### Variables
d_info = {}	# dict for entries in INFO field (all keys in the vcf get an entry)
format = []	# list for entries in FORMAT field
header=[]	# list for entries in header
isfirst = 1	# boolean isfirst - if true, we print the header once and turn off
number_samples = 0 	# number of samples in vcf file
# verboseheader = 1	# boolean verbose header

### Main
# read file from stdin
contents = sys.stdin.read().split("\n")

# loop through each line of the file
for line in contents: 

	# double # header 
	if ( line[0:2] == "##" ):
		# if argument, dont print double # header
		if (len(sys.argv) == 1): print(line)

		# e.g., double # header might look like this:
		# ##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
		# grab INFO subfield ID
		match = re.search(r'##INFO=<ID=(\w+)(.*)', line)
		# define d_info
		if match: 
			d_info[match.group(1)] = "-"

	# single # header
	elif ( line[0:1] == "#" ):
		# e.g., single # header looks like this:
		# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample_1        sample_2
		header = line.split("\t")
		# get the number of samples
		number_samples = len(header) - 9

	# the non-empty, non-header lines of the vcf
	elif ( line ):
		# if isfirst bool, print the header
		if (isfirst):
			# get list of elts in format field
			format = line.split("\t")[8].split(":")

			# print first part of header 
			print("\t".join(header[0:7]) + "\t"),

			# print INFO part of header -  sorted keys
			print("\t".join(sorted(d_info)) + "\t"),

			# print FORMAT part of header
			for i in range(1,1+number_samples): 
				myjoiner = "_" + str(i) + "\t"
				print(myjoiner.join(format) + "_" + str(i) + "\t"),

			print("\n"),

			# turn off flag
			isfirst = 0

		# reset dict
		for x in d_info: d_info[x]="-"

		# split line in tab
		linelist = line.split("\t")

		# print first part of line
		print("\t".join(linelist[0:7]) + "\t"),

		# loop thro INFO keys - if key found, add val
		for x in linelist[7].split(";"): 
			# pattern must be: blob=blog or we won't consider
			if (re.search(r'(\S+)=(\S+)', x)):
				d_info[x.split("=")[0]] = x.split("=")[1]
		# print INFO part of line - sorted values
		print("\t".join([d_info[x] for x in sorted(d_info)])),
		print("\t"),

		# print FORMAT part of header
		for i in range(9,9+number_samples): 
			print("\t".join(linelist[i].split(":")) + "\t"),

		print("\n"),
