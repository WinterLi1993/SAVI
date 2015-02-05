#!/usr/bin/env python

# About:
# This script takes a vcf file as input, piped in, and converts it to a tab-delimited text file, 
# with the aim to make the most human readable report

# Usage:
# cat myfile.vcf | vcf2report.py > myfile.txt
# To suppress lines beginning with "##" (i.e., the vcf header), run:
# cat myfile.vcf | vcf2report.py 0 > myfile.txt

import sys
import re

### Variables
d_info = {'EFF.EFFECT':'-'}	# dict for entries in INFO field (all keys in the vcf get an entry)
format = []			# list for entries in FORMAT field
header=[]			# list for entries in header
isfirst = 1			# boolean isfirst - if true, we print the header once and turn off
number_samples = 0 		# number of samples in vcf file
# verboseheader = 1		# boolean verbose header
num_eff_fields = 14		# number of SnpEff fields delimited by "|"
emptyefflist = ['-']*num_eff_fields	# define empty eff list 

### Main
# read file from stdin
contents = sys.stdin.read().split("\n")

# loop through each line of the file
for line in contents: 

	# double # header 
	if ( line[0:2] == "##" ):
		# if argument, print double # header
		if (len(sys.argv) == 1): print(line)

		# e.g., double # header might look like this:
		# ##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
		# grab INFO subfield ID
		match = re.search(r'##INFO=<ID=(\w+)(.*)', line)

		# define d_info
		if match: 
			# special case for EFF header, which has "|" delimited sub-fields
			if (match.group(1) == "EFF"):
				# eff line looks something like this:
				# ,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )' ">
				# get eff fields into header
				match2 = re.search(r'(.*)\((.*)\[ \| (.*)\] \)', line) 
				d_info.update(dict.fromkeys(["EFFECT%s.%s" % (chr(ord('a') + my_index), j.strip()) for my_index, j in enumerate( match2.group(2).split("|") + match2.group(3).split("|") )], "-"))
			else:
				d_info[match.group(1)] = "-"

	# single # header
	elif ( line[0:1] == "#" ):
		# e.g., single # header looks like this:
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		# but only want subset:
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		header = [v for j, v in enumerate(line.split("\t")) if j not in [5,6]]

		# get the number of samples (# cols after the FORMAT field)
		number_samples = len(header) - 7

	# the non-empty, non-header lines of the vcf
	elif ( line ):
		# if isfirst bool, print the header
		if (isfirst):
			# get list of elts in format field
			# but only subset: filter from
			# ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] to
			# ['SDP', 'RD', 'AD', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']
			format = [v for j, v in enumerate( line.split("\t")[8].split(":") ) if j not in [0,1,3,6,7]]

			# print first part of header, up until but not including INFO field
			print("\t".join(header[0:5]) + "\t"),

			# print INFO part of header -  sorted keys
			# print("\t".join(sorted(d_info)) + "\t"),
			
			# print selective parts of the INFO field
			for myfield in sorted(d_info): 
				# match savi fields, which can look like this: PD21_F, P1_F, P2_F, S1_P, S1_PF, S2_P, S2_PF, ...
				match4 = re.search(r'S(\d+)_(.*)', myfield) or re.search(r'P(\d+)_(.*)', myfield)

				if match4 or myfield.startswith('PD') or myfield == 'RS' or myfield == 'dbSNPBuildID' or myfield == 'COSMIC_NSAMP' or myfield == 'STRAND' or myfield == 'AA' or myfield == 'CDS' or myfield == 'NMutPerID' or myfield == 'MEGANORMAL_ID' or myfield == 'TOT_COUNTS' or myfield == 'SOM_VERIFIED_COUNTS' or myfield == 'SRC':
					print(myfield + "\t"),
				elif myfield.startswith('EFF'):
					print(myfield.split(".")[1] + "\t"),

			# print FORMAT part of header
			for i in range(1,1+number_samples): 
				myjoiner = "_" + str(i) + "\t"
				print(myjoiner.join(format) + "_" + str(i) + "\t"),

			print("\n"),

			# turn off flag
			isfirst = 0

		# reset dict
		for x in d_info: d_info[x]="-"

		# split line on tab
		linelist = line.split("\t")

		# now we want to delete QUAL FILTER from this:
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		del linelist[5:7]

		# print first part of line, up until but not including INFO field
		print("\t".join(linelist[0:5]) + "\t"),

		# loop thro INFO keys - if key found, add val
		for x in linelist[5].split(";"): 
			# pattern must be: blob=blob or we wont consider
			if (re.search(r'(\S+)=(\S+)', x)):
				d_info[x.split("=")[0]] = x.split("=")[1]
		# if eff key not present, add it
		if (not 'EFF' in d_info): 
			d_info['EFF'] = "-"

		## Keep these special keys in the INFO field:
		# COSMIC_NSAMP How many samples have this mutation
		# STRAND Gene strand
		# GENE Gene name
		# AA Peptide annotation
		# CDS CDS annotation
		# NMutPerID Number of mutations in the sample given by MEGANORMAL_ID
		# MEGANORMAL_ID disease|date|investigator|ID
		# IMPACT Impact of Mutation
		# TOT_COUNTS Total Mutation Count
		# SOM_VERIFIED_COUNTS Somatic Verified Mutation Count
		# GERM_COUNTS Germline Mutation Count
		# SOM_COUNTS Somatic Mutation Count
		# SRC Source List of Mutations

		# if key not present, add it
		if (not 'COSMIC_NSAMP' in d_info): 
			d_info['COSMIC_NSAMP'] = "-"
		if (not 'STRAND' in d_info): 
			d_info['STRAND'] = "-"
		if (not 'AA' in d_info): 
			d_info['AA'] = "-"
		if (not 'CDS' in d_info): 
			d_info['CDS'] = "-"
		if (not 'NMutPerID' in d_info): 
			d_info['NMutPerID'] = "-"
		if (not 'MEGANORMAL_ID' in d_info): 
			d_info['MEGANORMAL_ID'] = "-"
		if (not 'TOT_COUNTS' in d_info): 
			d_info['TOT_COUNTS'] = "-"
		if (not 'SOM_VERIFIED_COUNTS' in d_info): 
			d_info['SOM_VERIFIED_COUNTS'] = "-"
		if (not 'SRC' in d_info): 
			d_info['SRC'] = "-"
		if (not 'RS' in d_info): 
			d_info['RS'] = "-"
		if (not 'dbSNPBuildID' in d_info): 
			d_info['dbSNPBuildID'] = "-"

		# print INFO part of line - sorted values
		for myfield in sorted(d_info):
			# match savi fields
			match4 = re.search(r'S(\d+)_(.*)', myfield) or re.search(r'P(\d+)_(.*)', myfield)

			# SnpEff field is special case
			# It might look like:
			# INTRAGENIC(MODIFIER|||||AL583842.3||NON_CODING|||1)
			# SPLICE_SITE_REGION(LOW||||259|AL627309.1|protein_coding|CODING|ENST00000423372|1|1|WARNING_TRANSCRIPT_NO_START_CODON)
			if myfield == "EFF":
				myefflist = ['-']*num_eff_fields
				if (d_info[myfield] == "-"):
					print("\t".join(myefflist) + "\t"),
				else: 
					# split eff field on commas
					for ii in d_info[myfield].split(","):
						match3 = re.search(r'(\w+)\((.*)\)', ii)
						my_tmp_efflist = [match3.group(1)] + match3.group(2).split("|")
						# if subfield empty, set to '-'
						my_tmp_efflist = [myelt if myelt else '-' for myelt in my_tmp_efflist]
						# if no warnings or errors (list less than 13 elts), extend list to 13 elts
						mydiff = num_eff_fields - len(my_tmp_efflist)
						if (mydiff>0):
							my_tmp_efflist.extend(['-']*mydiff)
						# if efflist null, set it; else concat
						if (myefflist == emptyefflist):
							myefflist = my_tmp_efflist
						else:
							myefflist = ["%s,%s" % (jj, kk) for jj, kk in zip(myefflist, my_tmp_efflist)]
					# print output
					print("\t".join(myefflist) + "\t"),

			# else if match savi fields OR any special keys described above, print
			elif match4 or myfield.startswith('PD') or myfield == 'RS' or myfield == 'dbSNPBuildID' or myfield == 'COSMIC_NSAMP' or myfield == 'STRAND' or myfield == 'AA' or myfield == 'CDS' or myfield == 'NMutPerID' or myfield == 'MEGANORMAL_ID' or myfield == 'TOT_COUNTS' or myfield == 'SOM_VERIFIED_COUNTS' or myfield == 'SRC':
				print(d_info[myfield] + "\t"),
			
		# print FORMAT part of header
		for i in range(7,7+number_samples): 
			# but only want subset: filter from
			# ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] to
			# ['SDP', 'RD', 'AD', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']
			tmpformat = [v for j, v in enumerate(linelist[i].split(":")) if j not in [0,1,3,6,7]]
			print("\t".join(tmpformat) + "\t"),

		print("\n"),
