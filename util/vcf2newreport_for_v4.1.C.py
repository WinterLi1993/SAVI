#!/usr/bin/env python

### About:
# This script takes a vcf file as input, piped in, and converts it to a tab-delimited text file, 
# with the aim to make the most human readable report

### Usage:
# cat myfile.vcf | $0 > myfile.txt

import argparse
import sys
import re

### Variables
d_savi_info = {}			# dict for entries in INFO field pertaining to savi
d_samples = {}				# dict for sample names
format = []				# list for entries in FORMAT field
header=[]				# list for entries in header
isfirst = 1				# boolean isfirst - if true, we print the header once and turn off
number_samples = 0			# number of samples in vcf file
# verboseheader = 1			# boolean verbose header
# num_eff_fields = 14			# number of SnpEff fields delimited by "|"
num_eff_fields = 11			# number of SnpEff fields delimited by "|" (change to 11 to ignore GENOTYPE ERRORs and WARNINGS field)
emptyefflist = ['-']*num_eff_fields	# define empty eff list 

# make a human readable version:
header_readable = ['#chromosome', 'position', 'ref', 'alt', 'id.snp', 'snp_allele_freq', 'id.cosmic', 'cosmic_n.samples-mut']

# dict for entries in INFO field pertaining to SnpEff
d_eff_info = { 	'EFFECTa.Effect_Impact': '-',
		'EFFECTb.Functional_Class': '-',
		'EFFECTc.Codon_Change': '-',
		'EFFECTd.Amino_Acid_Change': '-',
		'EFFECTe.Amino_Acid_length': '-',
		'EFFECTf.Gene_Name': '-',
		'EFFECTg.Transcript_BioType': '-',
		'EFFECTh.Gene_Coding': '-',
		'EFFECTi.Transcript_ID': '-',
		'EFFECTj.Exon_Rank': '-' } 
		# 'EFFECTk.Genotype_Number': '-' } 

# dict for entries in INFO field - specify the INFO fields we want
d_info = { 	'Sgt1MAXFREQ': '-',
		'S1ADPP': '-',
		'COSMIC_NSAMP': '-',
		'TOT_COUNTS': '-',
		'SOM_VERIFIED_COUNTS': '-',
		'CAF': '-',
		'CNT': '-' } 

# dict to map vcf, SnpEff terms to human readable terms
d_readable = {	'#CHROM':'#chromosome', 
		'POS':'position',
		'ID':'id',
		'REF':'ref',
		'ALT':'alt',
		'SDP':'totdepth',
		'RD':'refdepth',
		'AD':'altdepth',
		'RBQ':'ref_ave_qual',
		'ABQ':'alt_ave_qual',
		'RDF':'ref_forward_depth',
		'RDR':'ref_reverse_depth',
		'ADF':'alt_forward_depth',
		'ADR':'alt_reverse_depth',
		'_P':'_presence_bool',
		'_PF':'_presence_posterior',
		'_F':'_freq',
		'_L':'_savi_change_stat',
		'_U':'_freq_upper',
		'COSMIC_NSAMP':'cosmic_n.samples-gene',
		'MEGANORMAL_ID':'meganormal_id',
		'NMutPerID':'n.mutations_per_meganormal_id',
		'RS':'reference_SNP_id(RS)',
		'SOM_VERIFIED_COUNTS':'cbio_somatic_verified_mutation_count',
		'SRC':'meganormal_186_TCGA_source',
		'TOT_COUNTS':'cbio_total_mutation_count',
		'dbSNPBuildID':'dbSNPBuild_for_RS',
		'STRAND':'strand',
		'CNT':'cosmic_n.samples-mut',
		'CAF':'snp_allele_freq',
		'Sgt1MAXFREQ':'Sgt1_max_frequency',
		'S1ADPP':'S1_alt_depth_per_position' }

# specify the format fields we want - only a subset of total:
# ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] to
format = ['SDP', 'RD', 'AD', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] 
# make a human readable version:
format_readable = [d_readable[j] for j in format]

# a function to test if INFO field is from SAVI - only want PD _L and P _F
def is_savi_info_field(mystring):
	"""Test if INFO field is from SAVI"""

	# grab savi specific INFO fields - e.g., 
	# no: ##INFO=<ID=PD21_F,Number=1,Type=Integer,Description="Savi freq delta for samples 2 vs 1">
	# yes: ##INFO=<ID=PD21_L,Number=1,Type=Integer,Description="Savi freq delta lower bound for samples 2 vs 1">
	# yes: ##INFO=<ID=PD21_U,Number=1,Type=Integer,Description="Savi freq delta upper bound for samples 2 vs 1">
	# yes: ##INFO=<ID=P1_F,Number=1,Type=Integer,Description="Savi freq for sample 1">
	# no: ##INFO=<ID=P1_L,Number=1,Type=Integer,Description="Savi freq lower bound for sample 1">
	# no: ##INFO=<ID=P1_U,Number=1,Type=Integer,Description="Savi freq upper bound for sample 1">
	# no: ##INFO=<ID=S1_P,Number=1,Type=Integer,Description="Savi presence call boolean for sample 1">
	# no: ##INFO=<ID=S1_PF,Number=1,Type=Float,Description="Savi posterior for presence or absence for sample 1">

	if re.search(r'P(\d)_F', mystring) or re.search(r'PD(\d+)_[UL]', mystring):
		return 1
	else:
		return 0

### Read Arguments
prog_description = """This script takes a vcf file as input, piped in, and converts it to a tab-delimited text file, 
with the aim to make a human readable report
"""

# http://docs.python.org/2/howto/argparse.html
parser = argparse.ArgumentParser(description=prog_description)
parser.add_argument("-v", "--verbose", 		action="store_true",		help="verbose mode. Default: off")
parser.add_argument("--header", 		action="store_true",		help="print vcf header lines beginning with \"##\". Default: off")
parser.add_argument("--debug", 			action="store_true",		help="turn on debugger (print extra stuff). Default: off")
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

		if match and is_savi_info_field(match.group(1)):
			d_savi_info[match.group(1)] = "-"

	# single # header
	elif ( line[0:1] == "#" ):
		# e.g., single # header looks like this:
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample_1', 'sample_2']

		header = line.split("\t")
		# get the number of samples (# cols after the FORMAT field)
		number_samples = len(header) - 9

		# check for error - did user give the correct number of sample names?
		if args.samples and len(d_samples) != number_samples:
			print("Error: Number of samples (" + str(number_samples) + ") doesn\'t match number of sample names the user provided (" + str(len(d_samples)) +")")
			sys.exit(0)

	# the non-empty, non-header lines of the vcf
	elif ( line ):
		# if isfirst bool - i.e., if we're on the first non-header line of the vcf, print the header of the report
		if (isfirst):

			if (args.debug):
				print("---DEBUG---")
				print(d_info)
				print(d_eff_info)
				print(d_savi_info)

			# print first part of header, up until but not including INFO field
			print("\t".join(header_readable) + "\t"),

			# print the SnpEff INFO part of header
			print("Effect\t"),
			for myfield in sorted(d_eff_info):

				print(myfield.split(".")[1] + "\t"),

			# print the non-savi INFO part of header
			for myfield in sorted(d_info):

				# make exceptions for CNT and CAF
				if myfield in d_readable and myfield != 'CNT' and myfield != 'CAF': 
					print(d_readable[myfield] + "\t"),

			# print the savi INFO part of header (selective parts of the INFO field)
			for myfield in sorted(d_savi_info):
				# match savi fields, which can look like this:
				# P1_F P2_F
				match4 = re.search(r'P(\d)(_F)', myfield)
				# or this:
				# PD21_L (assume fewer than 10 samples)
				match5 = re.search(r'PD(\d)(\d)(_L)', myfield)

				# make this field human readable
				if (match4): 
					# use sample names if user has provided them
					if (args.samples):
						print(d_samples[match4.group(1)] + d_readable[match4.group(2)] + "\t"),
					else:
						print(match4.group(1) + d_readable[match4.group(2)] + "\t"),

				# make this field human readable
				elif (match5): 
					# special case
					if ( match5.group(1) == "0" and match5.group(2) == "0" ):
						print("all-1" + d_readable[match5.group(3)] + "\t"),
					# use sample names if user has provided them
					elif (args.samples):
						print(d_samples[match5.group(1)] + "-" + d_samples[match5.group(2)] + d_readable[match5.group(3)] + "\t"),
					else:
						print(match5.group(1) + "-" + match5.group(2) + d_readable[match5.group(3)] + "\t"),

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

		# reset dicts
		for x in d_info: d_info[x]="-"
		for x in d_eff_info: d_eff_info[x]="-"
		for x in d_savi_info: d_savi_info[x]="-"

		# reset line string components
		# these variables will concatenate to form the complete line
		line_1 = ""
		line_2 = ""
		line_3 = ""
		line_4 = ""
		line_5 = ""
		line_6 = ""
		line_7 = ""
		line_8 = ""

		# split line on tab
		linelist = line.split("\t")

		# now we want to delete QUAL FILTER from this:
		#     0        1     2      3     4       5        6         7       8
		# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample_1', 'sample_2']
		del linelist[5:7]

		# print first part of line, not including ID, but otherwise up until but not including INFO field
		# print chrom pos
		line_1 = "\t".join(linelist[0:2]) + "\t"
		# print ref alt 
		line_2 = "\t".join(linelist[3:5]) + "\t"

		# get SNP portion of ID
		id_rs = [myid for myid in re.split('\W',linelist[2]) if myid.startswith('rs')]
		if id_rs:
			line_3 = ",".join(id_rs) + "\t"
		else:
			line_3 = "-\t"

		# get COSMIC portion of id
		id_cos = [myid for myid in re.split('\W',linelist[2]) if myid.startswith('C')]
		if id_cos:
			line_4 = ",".join(id_cos) + "\t"
		else:
			line_4 = "-\t"

		# start with SnpEff part empty by default
		line_5  = "\t".join(['-']*num_eff_fields) + "\t"

 		# loop thro INFO keys - if key found, add val
 		for x in linelist[5].split(";"): 
 			# pattern must be: blob=blob or we wont consider
 			if (re.search(r'(\S+)=(\S+)', x)):
 				# regular INFO field
 				if x.split("=")[0] in d_info:
 					d_info[x.split("=")[0]] = x.split("=")[1]

 				# SAVI INFO field
 				if is_savi_info_field(x.split("=")[0]):
 					d_savi_info[x.split("=")[0]] = x.split("=")[1]
 
 				# SnpEff INFO field is special case
 				# It might look like:
 				# INTRAGENIC(MODIFIER|||||AL583842.3||NON_CODING|||1)
 				# SPLICE_SITE_REGION(LOW||||259|AL627309.1|protein_coding|CODING|ENST00000423372|1|1|WARNING_TRANSCRIPT_NO_START_CODON)
 				if x.split("=")[0] == "EFF":

					# myefflist default empty
					myefflist = ['-']*num_eff_fields

					# split eff field on commas
					for ii in x.split("=")[1].split(","):
						match3 = re.search(r'(.*)\((.*)\)', ii)
						my_tmp_efflist = [match3.group(1)] + match3.group(2).split("|")
						# truncate to cut off errors and warnings
						my_tmp_efflist = my_tmp_efflist[:num_eff_fields]
						# if subfield empty, set to '-'
						my_tmp_efflist = [myelt if myelt else '-' for myelt in my_tmp_efflist]
						if (myefflist == emptyefflist):
							myefflist = my_tmp_efflist
						else:
							myefflist = ["%s,%s" % (jj, kk) for jj, kk in zip(myefflist, my_tmp_efflist)]
					# print output
					line_5 = "\t".join(myefflist) + "\t"

 		# print INFO part of line - sorted values
 		for myfield in sorted(d_info):
			if myfield != 'CNT' and myfield != 'CAF': 
 				line_6 += d_info[myfield] + "\t"

		# CNT and CAF special cases
		line_3 += d_info['CAF'] + "\t"
		line_4 += d_info['CNT'] + "\t"

 		# savi change string - this a composite of lower and upper bounds' vals 
 		savi_change = "nochange"
 		# print INFO part of line - sorted values
 		for myfield in sorted(d_savi_info):

			# this loop should look like:
			# P1_F P2_F P3_F PD21_L PD21_U PD32_L PD32_U

			if re.search(r'P(\d)_F', myfield):
 				line_7 += d_savi_info[myfield] + "\t"
			elif re.search(r'PD(\d+)_L', myfield):
				if (d_savi_info[myfield] != "-"):
					# lower bound > 0
					if int(d_savi_info[myfield]) > 0:
						savi_change = "up"
			elif re.search(r'PD(\d+)_U', myfield):
				if (d_savi_info[myfield] == "-"):
					savi_change = "-"
				# upper bound < 0
				elif int(d_savi_info[myfield]) < 0:
					savi_change = "down"
 				line_7 += savi_change + "\t"
 				# reset
 				savi_change = "nochange"

 		# print FORMAT part of header
 		for i in range(7,7+number_samples): 
 			# but only want subset: filter from
 			# ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR'] to
 			# ['SDP', 'RD', 'AD', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']
 			tmpformat = [v for j, v in enumerate(linelist[i].split(":")) if j not in [0,1,3,6,7]]
 			line_8 += "\t".join(tmpformat) + "\t"

		# print out line
		print(line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8)
