#!/usr/bin/env python

# About:
# This script takes a vcf file as input, piped in, and converts it to a tab-delimited text file, 
# with the aim to make the most human readable report

# Usage:
# cat myfile.vcf | $0 > myfile.txt

import argparse
import sys
import re

### Variables
d_info = {'EFF.EFFECT':'-'}		# dict for entries in INFO field (all keys in the vcf get an entry)
d_samples = {}				# dict for sample names
format = []				# list for entries in FORMAT field
header=[]				# list for entries in header
isfirst = 1				# boolean isfirst - if true, we print the header once and turn off
number_samples = 0			# number of samples in vcf file
# verboseheader = 1			# boolean verbose header
# num_eff_fields = 14			# number of SnpEff fields delimited by "|"
num_eff_fields = 12			# number of SnpEff fields delimited by "|" (change to 12 to ignore ERRORs and WARNINGS field)
emptyefflist = ['-']*num_eff_fields	# define empty eff list 

# dict to map vcf, SnpEff terms to human readable terms
d_readable = {'#CHROM':'#chromosome', 'POS':'position', 'ID':'id', 'REF':'ref', 'ALT':'alt', 'SDP':'totdepth', 'RD':'refdepth', 'AD':'altdepth', 'RBQ':'ref_ave_qual', 'ABQ':'alt_ave_qual', 'RDF':'ref_forward_depth', 'RDR':'ref_reverse_depth', 'ADF':'alt_forward_depth', 'ADR':'alt_reverse_depth', '_P':'_presence_bool', '_PF':'_presence_posterior', '_F':'_freq',  '_L':'_freq_lower_bool', '_U':'_freq_upper', 'COSMIC_NSAMP':'id.cosmic	cosmic_number_of_samples', 'MEGANORMAL_ID':'meganormal_id', 'NMutPerID':'number_of_mutations_per_meganormal_id', 'RS':'reference_SNP_id(RS)', 'SOM_VERIFIED_COUNTS':'cbio_somatic_verified_mutation_count', 'SRC':'meganormal_186_TCGA_source', 'TOT_COUNTS':'cbio_total_mutation_count', 'dbSNPBuildID':'dbSNPBuild_for_RS', 'STRAND':'strand', 'CAF':'id.snp	snp_allele_freq', 'Sgt1MAXFREQ':'Sgt1_max_frequency', 'S1ADPP':'S1_alt_depth_per_position'}

### Read Arguments
prog_description = """This script takes a vcf file as input, piped in, and converts it to a tab-delimited text file, 
with the aim to make a human readable report
"""

# http://docs.python.org/2/howto/argparse.html
parser = argparse.ArgumentParser(description=prog_description)
parser.add_argument("-v", "--verbose", 		action="store_true",		help="verbose mode. Default: off")
parser.add_argument("--header", 		action="store_true",		help="print vcf header lines beginning with \"##\". Default: off")
parser.add_argument("-s", "--samples",						help="comma delimited sample names")

args = parser.parse_args()

# create dict, d_samples, mapping indicies to user-provided sample names 
if (args.samples):
	# map indicies to sample names
	d_samples = {str(k+1) : j for k, j in enumerate(args.samples.split(','))}	

	# did user give more than 9 names? This is not yet supported
	if len(d_samples) > 9:
		print("Error: Number of samples greater than 9. The \"--samples\" flag does not support this. As an alternative, use the hack-y util/make_header.py")
		sys.exit(0)

### Main
# read file from stdin
contents = sys.stdin.read().split("\n")

# loop through each line of the file
for line in contents: 

	# if line is double # header, accumulate some dicts
	if ( line[0:2] == "##" ):
		# if argument, print double # header
		if (args.header): print(line)

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
				# d_info.update(dict.fromkeys(["EFFECT%s.%s" % (chr(ord('a') + my_index), j.strip()) for my_index, j in enumerate( match2.group(2).split("|") + match2.group(3).split("|") )], "-"))
				d_info.update(dict.fromkeys(["EFFECT%s.%s" % (chr(ord('a') + my_index), j.strip()) for my_index, j in enumerate( match2.group(2).split("|") )], "-"))
			else:
				d_info[match.group(1)] = "-"

	# single # header
	elif ( line[0:1] == "#" ):
		# e.g., single # header looks like this:
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		# but only want subset:
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		header = [v for j, v in enumerate(line.split("\t")) if j not in [5,6]]
		# make a human readable version:
		header_readable = ['#chromosome', 'position', 'ref', 'alt']

		# get the number of samples (# cols after the FORMAT field)
		number_samples = len(header) - 7

		# check for error - did user give the correct number of sample names?
		if args.samples and len(d_samples) != number_samples:
			print("Error: Number of samples (" + str(number_samples) + ") doesn\'t match number of sample names the user provided (" + str(len(d_samples)) +")")
			sys.exit(0)

	# the non-empty, non-header lines of the vcf
	elif ( line ):
		# if isfirst bool - i.e., if we're on the first non-header line of the vcf - finish printing the header of the report
		if (isfirst):
			# get list of elts in FORMAT field
			# but only subset: filter from
			# ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] to
			# ['SDP', 'RD', 'AD', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']
			format = [v for j, v in enumerate( line.split("\t")[8].split(":") ) if j not in [0,1,3,6,7]]
			# make a human readable version:
			format_readable = [d_readable[j] for j in format]

			# print first part of header, up until but not including INFO field
			print("\t".join(header_readable) + "\t"),

			# print INFO part of header (selective parts of the INFO field)
			for myfield in sorted(d_info): 
				# match savi fields, which can look like this:
				# P1_F P2_F
				match4 = re.search(r'([P])(\d)(_[FP](F?))', myfield)
				# or this:
				# PD21_L (assume fewer than 10 samples)
				match5 = re.search(r'(PD)(\d)(\d)(_[L])', myfield)

				# make this field human readable
				if (match4): 
					# use sample names if user has provided them
					if (args.samples):
						print(d_samples[match4.group(2)] + d_readable[match4.group(3)] + "\t"),
					else:
						print(match4.group(2) + d_readable[match4.group(3)] + "\t"),

				# make this field human readable
				elif (match5): 
					# use sample names if user has provided them
					if (args.samples):
						print(d_samples[match5.group(2)] + "-" + d_samples[match5.group(3)] + d_readable[match5.group(4)] + "\t"),
					else:
						print(match5.group(2) + "-" + match5.group(3) + d_readable[match5.group(4)] + "\t"),

				elif myfield == 'TOT_COUNTS' or myfield == 'SOM_VERIFIED_COUNTS' or myfield == 'CAF' or myfield == 'COSMIC_NSAMP' or myfield == 'Sgt1MAXFREQ' or myfield == 'S1ADPP':
					if myfield in d_readable:
						print(d_readable[myfield] + "\t"),
					else:
						print(myfield + "\t"),

				elif myfield.startswith('EFF'):
					print(myfield.split(".")[1] + "\t"),

			# print FORMAT part of header
			for i in range(1,1+number_samples): 
				# use sample names if user has provided them
				if (args.samples):
					myjoiner = "_" + d_samples[str(i)] + "\t"
					print(myjoiner.join(format_readable) + "_" + d_samples[str(i)] + "\t"),
				else:
					myjoiner = "_" + str(i) + "\t"
					print(myjoiner.join(format_readable) + "_" + str(i) + "\t"),

			print("\n"),

			# turn off flag
			isfirst = 0

		### now we're no longer on the first line:	

		# reset dict
		for x in d_info: d_info[x]="-"

		# split line on tab
		linelist = line.split("\t")

		# now we want to delete QUAL FILTER from this:
		#     0        1     2      3     4       5        6         7       8
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		del linelist[5:7]

		# print first part of line, not including ID, but otherwise up until but not including INFO field
		# print chrom pos
		print("\t".join(linelist[0:2]) + "\t"),
		# print ref alt 
		print("\t".join(linelist[3:5]) + "\t"),

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
		if (not 'Sgt1MAXFREQ' in d_info): 
			d_info['Sgt1MAXFREQ'] = "-"
		if (not 'S1ADPP' in d_info): 
			d_info['S1ADPP'] = "-"
		if (not 'COSMIC_NSAMP' in d_info): 
			d_info['COSMIC_NSAMP'] = "-"
		# if (not 'STRAND' in d_info): 
		#	d_info['STRAND'] = "-"
		# if (not 'AA' in d_info): 
		# 	d_info['AA'] = "-"
		# if (not 'CDS' in d_info): 
		# 	d_info['CDS'] = "-"
		# if (not 'NMutPerID' in d_info): 
		# 	d_info['NMutPerID'] = "-"
		# if (not 'MEGANORMAL_ID' in d_info): 
		# 	d_info['MEGANORMAL_ID'] = "-"
		if (not 'TOT_COUNTS' in d_info): 
			d_info['TOT_COUNTS'] = "-"
		if (not 'SOM_VERIFIED_COUNTS' in d_info): 
			d_info['SOM_VERIFIED_COUNTS'] = "-"
		if (not 'CAF' in d_info): 
			d_info['CAF'] = "-"
		# if (not 'SRC' in d_info): 
		# 	d_info['SRC'] = "-"
		# if (not 'RS' in d_info): 
		# 	d_info['RS'] = "-"
		# if (not 'dbSNPBuildID' in d_info): 
		# 	d_info['dbSNPBuildID'] = "-"

		# print INFO part of line - sorted values
		for myfield in sorted(d_info):
			# match savi fields
			match4 = re.search(r'P(\d+)_F', myfield)
			match5 = re.search(r'PD(\d+)_L', myfield)

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
						match3 = re.search(r'(.*)\((.*)\)', ii)
						my_tmp_efflist = [match3.group(1)] + match3.group(2).split("|")
						# truncate to cut off errors and warnings
						my_tmp_efflist = my_tmp_efflist[:num_eff_fields]
						# if subfield empty, set to '-'
						my_tmp_efflist = [myelt if myelt else '-' for myelt in my_tmp_efflist]
						# if no warnings or errors (list less than 13 elts), extend list to 13 elts
						# mydiff = num_eff_fields - len(my_tmp_efflist)
						# if (mydiff>0):
						# 	my_tmp_efflist.extend(['-']*mydiff)
						# if efflist null, set it; else concat
						if (myefflist == emptyefflist):
							myefflist = my_tmp_efflist
						else:
							myefflist = ["%s,%s" % (jj, kk) for jj, kk in zip(myefflist, my_tmp_efflist)]
					# print output
					print("\t".join(myefflist) + "\t"),

			# else if match savi fields OR any special keys described above, print
			elif match4 or myfield == 'TOT_COUNTS' or myfield == 'SOM_VERIFIED_COUNTS' or myfield == 'Sgt1MAXFREQ' or myfield == 'S1ADPP':
				print(d_info[myfield] + "\t"),
			elif myfield == 'CAF':
				# hack-ish: weld SNP column onto CAF column
				# get SNP portion of ID
				id_rs = [myid for myid in re.split('\W',linelist[2]) if myid.startswith('rs')]
				if id_rs:
					print(",".join(id_rs) + "\t"),
				else:
					print("-\t"),
				print(d_info[myfield] + "\t"),
			elif myfield == 'COSMIC_NSAMP':
				# hack-ish: weld COSMIC column onto COSMIC_NSAMP column
				# get COSMIC portion of id
				id_cos = [myid for myid in re.split('\W',linelist[2]) if myid.startswith('C')]
				if id_cos:
					print(",".join(id_cos) + "\t"),
				else:
					print("-\t"),
				print(d_info[myfield] + "\t"),
			elif match5:
				if (int(d_info[myfield]) < 0):
					print("-1\t"),
				elif (int(d_info[myfield]) > 0):
					print("1\t"),
				else:
					print("0\t"),
			
		# print FORMAT part of header
		for i in range(7,7+number_samples): 
			# but only want subset: filter from
			# ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] to
			# ['SDP', 'RD', 'AD', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']
			tmpformat = [v for j, v in enumerate(linelist[i].split(":")) if j not in [0,1,3,6,7]]
			print("\t".join(tmpformat) + "\t"),

		print("\n"),
