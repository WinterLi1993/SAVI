#!/usr/bin/env python

import argparse
import sys 
import re
import subprocess
import os
import itertools

# homemade:
from runsys.runsys import whichcmd, escape_special_char 

# -------------------------------------

# global variables
script_dir = os.path.dirname(__file__)
cwdir = os.getcwd()
# flag for make_qvt
qvtargs = ""
# prior_dict = {}

def main():
	"""Main block"""

	args = get_arg()
	savidir, inplist = run_setup(args)
	run_vlad_code_present(args, savidir, inplist)
	run_vlad_code_freq(args, savidir, inplist)
	run_vlad_code_compare(args, savidir, inplist)

def compress_vcf(args, myvcf):
	"""Compress a vcf file with bgzip and tabix"""

	# make sure a zipped version of the file doesnt already exist
	cmd = "rm -f " + myvcf + ".bgz"
	whichcmd(cmd, args, 0)

	cmd = "bgzip " + myvcf
	whichcmd(cmd, args, 0)

	# move extension from .gz to .bgz
	cmd = "mv " + myvcf + ".gz " + myvcf + ".bgz "
	whichcmd(cmd, args, 0)

	cmd = "tabix -p vcf " + myvcf + ".bgz"
	whichcmd(cmd, args, 0)
	
def run_setup(args):
	"""Run preliminary set up"""

	# make dirs
	cmd="mkdir -p " + args.outputdir
	if (args.verbose): print("\n# make directories");
	whichcmd(cmd, args, 0)
	
	savidir = args.outputdir
	
	# input list is the list of uniq sample indices (split on either "," or ":")
	inplist = list(set(re.split(",|:", args.sample)))
	inplist.sort()

	if ( len(prior_dict.keys()) != len(inplist) ):
		print("Please supply a prior for every sample in your list:")
		print("input: "),
		print(inplist)
		print("prior: "),
		print(prior_dict.keys())
		sys.exit(1)

	return savidir, inplist

def run_vlad_code_present(args, savidir, inplist):
	"""Run Vladimir's savi computations to check for presence and filter based on it"""

	# declare as global or else Python will treat as a local variable
	global qvtargs

	# variable to hold comma-delimited list of SGE job IDs
	sgejobids=""

	# this variable is true only for the first iteration of the list
	bool_first = 1
	# a string to be used by vcffilter - we'll build it up as we go
	filter_str = ""

	# prepend different things to the command string, depending on if whole file or small region
	prependcmd = "" 

	if (not args.region):
		prependcmd = "zcat " + args.input + " | "
	else:
		prependcmd = "tabix -h " + args.input + " " + args.region + " | "

	# make sure the header file doesnt exist because we're going to append to it
	cmd = "rm -f " + savidir + "/{header_addition.txt,vcfheader.txt}"
	whichcmd(cmd, args, 0)

	# compute presence on everything in the input list (this shell scripting is a little messy)
	for i in inplist:
		if bool_first:
			# create a string to give to vcffilter ( vcffilter doesnt support less than or equals to <= )
			filter_str = "( S" + i + "_P = 1 & " + "S" + i + "_PF < " + args.savipresent + " )"
			bool_first = 0

			# if hybrid, start with tot depth = RD + AD for first iteration, then turn off
			if ( args.hybrid ): qvtargs = "--rdplusad"
	
		else:
			# update string for sample in inplist
			filter_str = filter_str + " | ( S" + i + "_P = 1 & " + "S" + i + "_PF < " + args.savipresent + " )"

			# if hybrid, start with tot depth = RD + AD for first iteration, then turn off
			if ( args.hybrid ): qvtargs = ""

		# make conf files

		# cat testout/sample_1_1.1p | ../savi/savi_poster -p testout/prior | ../savi/savi_conf -z | cut -f4- > testout/savi/conf_1.txt
		# cat testout/sample_2_1.1p | ../savi/savi_poster -p testout/prior | ../savi/savi_conf -z | cut -f4- > testout/savi/conf_2.txt
		# cat testout/sample_3_1.1p | ../savi/savi_poster -p testout/prior | ../savi/savi_conf -z | cut -f4- > testout/savi/conf_3.txt
		# cat testout/sample_4_1.1p | ../savi/savi_poster -p testout/prior | ../savi/savi_conf -z | cut -f4- > testout/savi/conf_4.txt

		if (args.verbose): print("\n# make a presence call for sample " + i);
		cmd = prependcmd + args.bin + "/make_qvt -1 -s " + i  + " " + qvtargs + " | " + \
		args.bin + "/savi_poster -p " + prior_dict[i] + " | " + \
		args.bin + "/savi_conf -s " + args.saviprecision + " | cut -f4- | awk -v samp=" + i + " '" + '{print "S"samp"_P="$1";S"samp"_PF="$2}' + "' > " + savidir + "/conf_" + i + ".txt"
		# args.bin + "/savi_conf -z | cut -f4- | awk -v samp=" + i + " '" + '{print "S"samp"_P="$1";S"samp"_PF="$2}' + "' > " + savidir + "/conf_" + i + ".txt"
		# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
		# run cmd and store SGE job id
		cmd=escape_special_char(cmd)
		myjobid = whichcmd(cmd, args, 1, args.sge, "j_s" + i + "_pres")

		if ( args.sge ):
			print("Your job " + myjobid + " has been submitted")
			sgejobids = myjobid + "," + sgejobids

		# fix vcf header
		with open(savidir + "/header_addition.txt", 'a') as f:
			f.write("##INFO=<ID=S" + i + "_P,Number=1,Type=Integer,Description=\"Savi presence call boolean for sample " + i + "\">\n")
			# scientific notation supported? Type=Float? this seems to work
			f.write("##INFO=<ID=S" + i + "_PF,Number=1,Type=Float,Description=\"Savi posterior for presence or absence for sample " + i + "\">\n")
	
	if (args.verbose): print("\n# paste savi presence numbers into INFO field");
	cmd = 'paste -d";" ' + savidir + "/conf_*.txt > " + savidir + "/conf.txt"
	cmd=escape_special_char(cmd)
	# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
	mysyncjob = whichcmd(cmd, args, 1, args.sge, "mysync", sgejobids, 1)

	if ( args.sge ):
		print("Your job " + mysyncjob + " has been submitted")

	if (not args.debug): 
		if (args.verbose): print("\n# clean up");
		cmd = "rm -f " + savidir + "/conf_*.txt"
		whichcmd(cmd, args, 0)

	# add into INFO field
	cmd = prependcmd + args.bin + "/add_to_info -f " + savidir + "/conf.txt --header " + savidir + "/vcfheader.txt > " + savidir + "/tmp0.txt"
	whichcmd(cmd, args, 0)

	# add new header lines to vcf
	cmd = "cat " + savidir + "/header_addition.txt " + savidir + "/vcfheader.txt " + savidir + "/tmp0.txt > " + savidir + "/addsavi.vcf"
	whichcmd(cmd, args, 0)

	if (not args.debug): 
		if (args.verbose): print("\n# clean up");
		cmd = "rm -f " + savidir + "/{header_addition.txt,vcfheader.txt,tmp0.txt,conf.txt}"
		whichcmd(cmd, args, 0)

	# don't filter for presence
	if (args.nofilter):
		# compress
		if (args.verbose): print("\n# compress addsavi.vcf");
		compress_vcf(args, savidir + "/addsavi.vcf")
	
	# filter for presence
	else:
		# filter with vcffilter 
		if (args.verbose): print("\n# filter for present variants");
		# does this work on bgzipped files?
		cmd = "vcffilter -f \"" + filter_str + "\" " + savidir + "/addsavi.vcf > " + savidir + "/filtersavi.vcf"
		cmd = escape_special_char(cmd)
		whichcmd(cmd, args, 0)

		if (not args.debug): 
			if (args.verbose): print("\n# clean up");
			cmd = "rm -f " + savidir + "/addsavi.vcf"
			whichcmd(cmd, args, 0)

		# compress
		if (args.verbose): print("\n# compress filtersavi.vcf");
		compress_vcf(args, savidir + "/filtersavi.vcf")

def run_vlad_code_freq(args, savidir, inplist):
	"""Run Vladimir's savi computations to get freq"""

	# declare as global or else Python will treat as a local variable
	global qvtargs

	# variable to hold comma-delimited list of SGE job IDs
	sgejobids=""

	# this variable is true only for the first iteration of the list
	bool_first = 1

	# make sure the header file doesnt exist because we're going to append to it
	cmd = "rm -f " + savidir + "/{header_addition.txt,vcfheader.txt}"
	whichcmd(cmd, args, 0)

	# prepend this to cmd str:
	firstcmd=""
	if (args.nofilter):
		firstcmd = "zcat " + savidir + "/addsavi.vcf.bgz | " 
	else:
		firstcmd = "zcat " + savidir + "/filtersavi.vcf.bgz | " 

	# make freq files
	for i in inplist:
		if bool_first:
			# if hybrid, start with tot depth = RD + AD for first iteration, then turn off
			if ( args.hybrid ): qvtargs = "--rdplusad"
			bool_first = 0
		else:
			# if hybrid, start with tot depth = RD + AD for first iteration, then turn off
			if ( args.hybrid ): qvtargs = ""

		if (args.verbose): print("\n# compute freq for sample " + i);
		cmd = firstcmd + args.bin + "/make_qvt -1 -s " + i + " " + qvtargs + " | " + \
		args.bin + "/savi_poster -p " + prior_dict[i] + " | " + \
		args.bin + "/savi_conf -fs " + args.saviconf + " " + args.saviprecision + " | awk -v samp=" + i + " '" + '{mystr="P"samp; print mystr"_F="$1";"mystr"_L="$5";"mystr"_U="$(NF-1)}' + "' > " + savidir + "/freq_" + i + ".txt"
		# args.bin + "/savi_conf -fc " + args.saviconf + " | awk -v samp=" + i + " '" + '{mystr="P"samp; print mystr"_F="$1";"mystr"_L="$5";"mystr"_U="$(NF-1)}' + "' > " + savidir + "/freq_" + i + ".txt"
		cmd = escape_special_char(cmd)
		# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
		# run cmd and store SGE job id
		myjobid = whichcmd(cmd, args, 1, args.sge, "j_s" + i + "_freq")

		if ( args.sge ):
			print("Your job " + myjobid + " has been submitted")
			sgejobids = myjobid + "," + sgejobids

		# fix vcf header
		with open(savidir + "/header_addition.txt", 'a') as f:
			f.write("##INFO=<ID=P" + i + "_F,Number=1,Type=Integer,Description=\"Savi freq for sample " + i + "\">\n")
			f.write("##INFO=<ID=P" + i + "_L,Number=1,Type=Integer,Description=\"Savi freq lower bound for sample " + i + "\">\n")
			f.write("##INFO=<ID=P" + i + "_U,Number=1,Type=Integer,Description=\"Savi freq upper bound for sample " + i + "\">\n")
	
	if (args.verbose): print("\n# paste savi freq numbers into INFO field");
	cmd = 'paste -d";" ' + savidir + "/freq_*.txt > " + savidir + "/freq.txt"
	cmd = escape_special_char(cmd)
	# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
	mysyncjob = whichcmd(cmd, args, 1, args.sge, "mysync", sgejobids, 1)

	if ( args.sge ):
		print("Your job " + mysyncjob + " has been submitted")

	if (not args.debug): 
		if (args.verbose): print("\n# clean up");
		cmd = "rm -f " + savidir + "/freq_*.txt"
		whichcmd(cmd, args, 0)

	# add into INFO field
	cmd = firstcmd + args.bin + "/add_to_info -f " + savidir + "/freq.txt --header " + savidir + "/vcfheader.txt > " + savidir + "/tmp0.txt"
	whichcmd(cmd, args, 0)

	# add new header lines to vcf
	cmd = "cat " + savidir + "/header_addition.txt " + savidir + "/vcfheader.txt " + savidir + "/tmp0.txt > " + savidir + "/freqsavi.vcf"
	whichcmd(cmd, args, 0)

	if (not args.debug): 
		if (args.verbose): print("\n# clean up");
		cmd = "rm -f " + savidir + "/{header_addition.txt,vcfheader.txt,tmp0.txt,freq.txt}"
		whichcmd(cmd, args, 0)

	# compress
	if (args.verbose): print("\n# compress freqsavi.vcf");
	compress_vcf(args, savidir + "/freqsavi.vcf")

def run_vlad_code_compare(args, savidir, inplist):
	"""Run Vladimir's savi computations to compare samples"""

	# declare as global or else Python will treat as a local variable
	global qvtargs

	# variable to hold comma-delimited list of SGE job IDs
	sgejobids=""

	# if hybrid, SAMP 1 tot depth = RD + AD, SAMP 2 tot depth = SDP
	if ( args.hybrid ): qvtargs = "--hybrid"

	# make sure the header file doesnt exist because we're going to append to it
	cmd = "rm -f " + savidir + "/{header_addition.txt,vcfheader.txt}"
	whichcmd(cmd, args, 0)

	if (args.verbose): print("\n# savi comparison");
	pairwise_list = [s for s in args.sample.split(",") if ":" in s]

	for s in pairwise_list:
		# the priors
		str_a = s.split(":")[0]
		str_b = s.split(":")[1]

		prior_a = prior_dict[str_a]
		prior_b = prior_dict[str_b]

		# paste testout/sample_1_1.1p testout/sample_2_1.1p | cut -f1-4,6-8 | ../savi/savi_poster -pd testout/prior testout/prior | ../savi/savi_conf -fc 1e-5  > testout/savi/pd_12.txt
		# paste testout/sample_1_1.1p testout/sample_3_1.1p | cut -f1-4,6-8 | ../savi/savi_poster -pd testout/prior testout/prior | ../savi/savi_conf -fc 1e-5  > testout/savi/pd_13.txt
		# paste testout/sample_1_1.1p testout/sample_4_1.1p | cut -f1-4,6-8 | ../savi/savi_poster -pd testout/prior testout/prior | ../savi/savi_conf -fc 1e-5  > testout/savi/pd_14.txt
		# paste testout/sample_2_1.1p testout/sample_3_1.1p | cut -f1-4,6-8 | ../savi/savi_poster -pd testout/prior testout/prior | ../savi/savi_conf -fc 1e-5  > testout/savi/pd_23.txt
		# paste testout/sample_2_1.1p testout/sample_4_1.1p | cut -f1-4,6-8 | ../savi/savi_poster -pd testout/prior testout/prior | ../savi/savi_conf -fc 1e-5  > testout/savi/pd_24.txt
		# paste testout/sample_3_1.1p testout/sample_4_1.1p | cut -f1-4,6-8 | ../savi/savi_poster -pd testout/prior testout/prior | ../savi/savi_conf -fc 1e-5  > testout/savi/pd_34.txt

		cmd = "zcat " + savidir + "/freqsavi.vcf.bgz | " + \
		args.bin + "/make_qvt -1 -2s " + str_a + "," + str_b + " " + qvtargs + " | " + \
		args.bin + "/savi_poster -pd " + prior_a + " " + prior_b + " | " + \
		args.bin + "/savi_conf -fs " + args.saviconf + " " + args.saviprecision + " | awk -v samp1=" + str_a + " -v samp2=" + str_b + " '" + '{mystr="PD"samp1 samp2; print mystr"_F="$1";"mystr"_L="$5";"mystr"_U="$(NF-1)}' + "' > " + savidir + "/pd_" + str_a + str_b + ".txt"
		# args.bin + "/savi_conf -fc " + args.saviconf + " | awk -v samp1=" + str_a + " -v samp2=" + str_b + " '" + '{mystr="PD"samp1 samp2; print mystr"_F="$1";"mystr"_L="$5";"mystr"_U="$(NF-1)}' + "' > " + savidir + "/pd_" + str_a + str_b + ".txt"
		cmd = escape_special_char(cmd)

		# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
		# run cmd and store SGE job id
		myjobid = whichcmd(cmd, args, 1, args.sge, "j_s" + str_a + "_s" + str_b + "_cmp")

		if ( args.sge ):
			print("Your job " + myjobid + " has been submitted")
			sgejobids = myjobid + "," + sgejobids

		# fix vcf header
		with open(savidir + "/header_addition.txt", 'a') as f:
			f.write("##INFO=<ID=PD" + str_a + str_b + "_F,Number=1,Type=Integer,Description=\"Savi freq delta for sample " + str_a + " vs " + str_b + "\">\n")
			f.write("##INFO=<ID=PD" + str_a + str_b + "_L,Number=1,Type=Integer,Description=\"Savi freq delta lower bound for sample " + str_a + " vs " + str_b + "\">\n")
			f.write("##INFO=<ID=PD" + str_a + str_b + "_U,Number=1,Type=Integer,Description=\"Savi freq delta upper bound for sample " + str_a + " vs " + str_b + "\">\n")

	# a there was at least one savi comparison and there are more than 2 samples, run a 1 vs ALL comparision per order of JiGuang
	if ( pairwise_list and len(inplist) > 2 ):
		# arbitarily use the first prior	
		str_a = pairwise_list[0].split(":")[0]
		str_b = pairwise_list[0].split(":")[0]

		prior_a = prior_dict[str_a]
		prior_b = prior_dict[str_b]

		cmd = "zcat " + savidir + "/freqsavi.vcf.bgz | " + \
		args.bin + "/make_qvt -1 -1vsall | " + \
		args.bin + "/savi_poster -pd " + prior_a + " " + prior_b + " | " + \
		args.bin + "/savi_conf -fs " + args.saviconf + " " + args.saviprecision + " | awk -v samp1=0 -v samp2=0" + " '" + '{mystr="PD"samp1 samp2; print mystr"_F="$1";"mystr"_L="$5";"mystr"_U="$(NF-1)}' + "' > " + savidir + "/pd_00.txt"
		cmd = escape_special_char(cmd)

		# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
		# run cmd and store SGE job id
		myjobid = whichcmd(cmd, args, 1, args.sge, "j_s" + str_a + "_s" + str_b + "_cmp")

		if ( args.sge ):
			print("Your job " + myjobid + " has been submitted")
			sgejobids = myjobid + "," + sgejobids

		# fix vcf header
		with open(savidir + "/header_addition.txt", 'a') as f:
			f.write("##INFO=<ID=PD00_F,Number=1,Type=Integer,Description=\"Savi freq delta for samples 1 vs all the others\">\n")
			f.write("##INFO=<ID=PD00_L,Number=1,Type=Integer,Description=\"Savi freq delta lower bound for samples 1 vs all the others\">\n")
			f.write("##INFO=<ID=PD00_U,Number=1,Type=Integer,Description=\"Savi freq delta upper bound for samples 1 vs all the others\">\n")

	# a there was at least one savi comparison
	if (pairwise_list):

		if (args.verbose): print("\n# paste savi numbers into INFO field");
		cmd = 'paste -d";" ' + savidir + "/pd_*.txt > " + savidir + "/pd.txt"
		cmd = escape_special_char(cmd)
		# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
		mysyncjob = whichcmd(cmd, args, 1, args.sge, "mysync", sgejobids, 1)

		if ( args.sge ):
			print("Your job " + mysyncjob + " has been submitted")

		if (not args.debug): 
			if (args.verbose): print("\n# clean up");
			cmd = "rm -f " + savidir + "/pd_*.txt"
			whichcmd(cmd, args, 0)

		cmd = "zcat " + savidir + "/freqsavi.vcf.bgz | " + \
		args.bin + "/add_to_info -f " + savidir + "/pd.txt --header " + savidir + "/vcfheader.txt > " + savidir + "/tmp0.txt"
		whichcmd(cmd, args, 0)

		cmd = "cat " + savidir + "/header_addition.txt " + savidir + "/vcfheader.txt " + savidir + "/tmp0.txt > " + savidir + "/finalsavi.vcf"
		whichcmd(cmd, args, 0)

		if (args.verbose): print("\n# compress finalsavi.vcf");
		compress_vcf(args, savidir + "/finalsavi.vcf")

		if (not args.debug): 
			if (not args.keepfreqfile):
				if (args.verbose): print("\n# clean up");
				cmd = "rm -f " + savidir + "/{freqsavi.vcf.bgz,freqsavi.vcf.bgz.tbi}"
				whichcmd(cmd, args, 0)

	if (not args.debug): 
		if (args.verbose): print("\n# clean up");
		cmd = "rm -f " + savidir + "/{header_addition.txt,vcfheader.txt,tmp0.txt,pd.txt,addsavi.vcf.bgz,addsavi.vcf.bgz.tbi,filtersavi.vcf.bgz,filtersavi.vcf.bgz.tbi}"
		whichcmd(cmd, args, 0)

	if ( args.verbose ): print("[END]")

def get_arg():
	"""Get Arguments"""

	prog_description = """Run savi pipeline. SAVI STEP.
	You must supply a (multi-sample) vcf and prior files for each sample in the vcf you want to run savi on.
	This script will take your vcf and filter for presence to produce a subset vcf.
	The script will then add the savi numbers to the INFO field of this filtered vcf.
	Notes: This wrapper needs to find the savi binaries in order to function properly. 
	It assumes they are in a directory called savi/ in the same folder where this script resides.
	If this is not the case, you can supply the path.
	You must also supply a directory for the script to work in. (If you don't, it will make a folder tmp/ in the cwd and work there)
	This script has the following dependencies: tabix and bgzip (available here: http://samtools.sourceforge.net/tabix.shtml) and vcflib (available here: https://github.com/ekg/vcflib)
	"""

	# http://docs.python.org/2/howto/argparse.html
	parser = argparse.ArgumentParser(description=prog_description)
	parser.add_argument("-v", "--verbose", 						action="store_true",		help="verbose mode. Default: off")
	parser.add_argument("--version", 						action="store_true",		help="print version and exit")
	parser.add_argument("--example", 						action="store_true",		help="print example command and exit")
	parser.add_argument("--debug", 							action="store_true",		help="dont clean up intermediate files. Default: off")
	# argument as well as pipe-able
	parser.add_argument("-i", "--input",										help="a (multi-sample) vcf file whose sample fields at least have the subfields SDP RBQ ABQ RD AD. Must be bgzipped and indexed for tabix")
	parser.add_argument("-r", "--region",										help="the region in your vcf you want to run savi on")
	parser.add_argument("-n", "--name",						default="sample",		help="the name of your sample. Default: sample")
	parser.add_argument("-s", "--sample",										help="a comma-delimited list of the sample indices in the vcf on which to run savi (use colons to denote pairwise comparisons (e.g., 1,3:2,4:2 would be 1, 3 minus 2, and 4 minus 2)). If you leave this blank, the script will make all pairwise comparisons")
	parser.add_argument("-p", "--prior",										help="a comma-delimited list of sample:prior for each sample in the vcf file you are running variants on (e.g., 1:file1,3:file3 if you are running savi on samples 1 and 3)")
	parser.add_argument("-b", "--bin",						default = script_dir + "/bin",	help="directory where the savi binaries reside. Default: [scriptdir]/bin")
	parser.add_argument("-o", "--outputdir",					default = cwdir + "/tmp", 	help="output directory. Default: [cwd]/tmp")
	# TO DO give these better names, reflecting what they actually are
	parser.add_argument("--nofilter", 						action="store_true",		help="do not filter according to savi presence. Default: off (i.e., filter)")
	parser.add_argument("--rdplusad", 						action="store_true",		help="Use AD + RD as tot depth instead of SDP. Default: off")
	parser.add_argument("--hybrid", 						action="store_true",		help="Use AD + RD as tot depth for sample 1, SDP as total depth for sample 2. Default: off")
	parser.add_argument("--saviconf",						default = "1e-5", 		help="the savi conf parameter. Default: 1e-5")
	parser.add_argument("--savipresent",						default = "1e-6", 		help="the savi presence parameter. Default: 1e-6")
	parser.add_argument("--saviprecision",						default = "0", 			help="the savi precision parameter added by Hossein. Default: 0")
	parser.add_argument("--keepfreqfile",						default = "0", 			help="keep the freq file. Default: 0")
	parser.add_argument("-q", "--sge",						action="store_true",		help="qsub commands with SGE. Default: off")
	parser.add_argument("-l", "--sgelog",						default = cwdir + "/log", 	help="for SGE only. Default: [cwd]/log")
	parser.add_argument("-m", "--sgemem",						default ="4", 			help="for SGE only. sge memory allotment (hrs). Default: 4")
	parser.add_argument("-t", "--sgetime",						default ="4", 			help="for SGE only. sge time allotment (Gigabytes). Default: 4")
	args = parser.parse_args()

	# check for any input
	if ( not len(sys.argv) > 1 ):
		parser.print_help()
		sys.exit(0)

	if ( args.version ):
		print( os.path.basename(__file__) + ": Version 1.0" )
		sys.exit(0)

	# print example command 
	if ( args.example ):
		my_example = " -v -i test.vcf.bgz -r 1:17000-18000 -n test -s 1,3:2,4:3 -p 1:dir_prior/prior1,2:dir_prior/prior2,3:dir_prior/prior3,4:dir_prior/prior4 -o . --savipresent 5e-01\n\n" + \
			     "This would compute savi numbers for rows in test.vcf.bgz in the region chr1:17000-18000 for sample 1, sample 3 minus sample 2, and sample 4 minus sample 3.\n" + \
			     "given a prior file for sample 1, dir_prior/prior1, a prior file for sample 2, dir_prior/prior2, a prior file for sample 3, dir_prior/prior3, and a \n" + \
			     "prior file for sample 4, dir_prior/prior4.\n" + \
			     "The output would be saved in the cwd"
		print( "# " + os.path.basename(__file__) + my_example )
		sys.exit(0)

	# check for input vcf
	if ( not args.input or not os.path.isfile(args.input) ):
		print("can't find input vcf")
		sys.exit(1)

	# check for input prior
	if ( not args.prior ):
		print("can't find prior")
		sys.exit(1)

	# check for priors
	#print("prior: "),
	#print(args.prior)
	mytmp=re.split(':|,', args.prior)
	global prior_dict
	prior_dict = dict(mytmp[i:i+2] for i in range(0, len(mytmp), 2))
	#print(prior_dict)
	for i in prior_dict.values():
	 	if ( not os.path.isfile(os.path.expanduser(i)) ):
	 		print("can't find prior " + i)
	 		sys.exit(1)

	# check for binaries
	for j in ("savi_poster", "savi_poster_accum", "savi_poster_merge", "make_qvt"):
		if ( not os.path.isfile(args.bin + "/" + j) ):
			print("can't find " + args.bin + "/" + j)
			sys.exit(1)

	# set qvtargs
	global qvtargs
	if ( args.rdplusad ):
		qvtargs = "--rdplusad"

	if ( args.verbose ):
		print("[BEGIN]")
		print("Arguments: "),
		print(args)

	return args

# -------------------------------------

if __name__ == "__main__":

	main()	
