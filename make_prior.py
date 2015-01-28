#!/usr/bin/env python

import argparse
import sys 
import re
import subprocess
import os
# import pylab

# homemade:
from runsys.runsys import whichcmd
from runsys.runsys import run_qsub

# -------------------------------------

script_dir = os.path.dirname(__file__)
cwdir = os.getcwd()

def main():
	"""Main block"""

	args = get_arg()
	run_vlad_code(args)

def run_setup(args):
	"""Run preliminary set up"""

	# make dirs
	cmd="mkdir -p " + args.outputdir + "/{inps,iters,res}"
	if (args.verbose): print("# make directories");
	whichcmd(cmd, args, 0)

	# make SGE dirs
	if ( args.sge ):
		cmd="mkdir -p " + args.sgelog
		whichcmd(cmd, args, 0)

	# location of dirs
	resdir = args.outputdir + "/res"
	itdir = args.outputdir + "/iters"
	inpdir = args.outputdir + "/inps"

	# number of regions over which to parallelize
	numregions = 1

	# copy prior to res dir
	cmd="cp " + args.prior + " " + resdir + "/1"
	if (args.verbose): print("# copy first prior to res/ dir"); 
	whichcmd(cmd, args, 0)

	if (args.verbose): print("# number of iterations: " + str(args.iteration));

	# if vcf
	if (args.input):
		if (args.region):
			numregions = len(args.region.split(","))
		else:
			if (args.verbose): print("# you did not specify any regions"); 
			# set args.region to empty
			args.region = ""
	# if qvt input piped in via stdin, make tmp 1p file
	elif (args.qvt):
		if (args.qvt == "-"):
			# location of inp file
			inpfile = args.outputdir + '/inps/1'
			cmd="rm -rf " + inpfile
			whichcmd(cmd, args, 0)

			if (args.verbose): print("# save copy of inp: " + inpfile); 
			with open(inpfile, 'w') as f:
				for line in sys.stdin:
					f.write(line)

			if ( os.path.getsize("myfile") > 0 ):
				print("can't find input")
				sys.exit(1)
		# else link qvt
		else:
			numregions = len(args.qvt.split(","))

			# rm preexisting input files
			cmd="rm -rf " +  " ".join([inpdir + "/%d" % x for x in range(1,numregions+1)])
			whichcmd(cmd, args, 0)

			if (args.verbose): print("# link inp file"); 
			for i, myinp in enumerate(args.qvt.split(",")):
				cmd="ln -sf " + myinp + " " + args.outputdir + '/inps/' + str(i+1)
				whichcmd(cmd, args, 0)

	return resdir, itdir, inpdir, numregions

def run_vlad_code(args):
	"""Run Vladimir's prior building binaries"""

	resdir, itdir, inpdir, numregions =  run_setup(args)

	# i is iterations of prior building
	# j is the number of (q,v,t) files parallelized over

	# variable to hold comma-delimited list of SGE job IDs
	sgejobids=""

	for i in range(1,args.iteration+1):

		if (args.verbose): print("iteration: " + str(i))

		# if after the first iteration
		if (i>1):
			mergefiles = " ".join([itdir + "/%d" % x for x in range(1,numregions+1)])
			cmd = args.bin + "/savi_poster_merge " + mergefiles + " > " + resdir + "/" + str(i)
			whichcmd(cmd, args, 0)
			cmd = "rm " + mergefiles
			# whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
			mysyncjob = whichcmd(cmd, args, 1, args.sge, "j_s" + str(args.sampleindex) + "_i" + str(i) + "_merge", sgejobids, 1)

			if ( args.sge ):
				print("Your job " + mysyncjob + " has been submitted")

			# clear sgejobids
			sgejobids=""

		# if less than the last iteration
		if (i<args.iteration):
			for j in range(1,numregions+1):
				# prepend to cmd
				prependcmd = ""
				if (args.input):
					if (args.region == ""):
						prependcmd = "zcat " + args.input + " | " + args.bin + "/make_qvt -1 -s " + str(args.sampleindex)  + " | "
					else:
						prependcmd = "tabix " + args.input + " " + args.region.split(",")[j-1] + " | " + args.bin + "/make_qvt -1 -s " + str(args.sampleindex)  + " | "
				elif (args.qvt):
					prependcmd = "cat " + args.outputdir + "/inps/" + str(j) + " | "

				cmd = prependcmd + args.bin + "/savi_poster -p " + resdir + "/" + str(i) + " | " + \
				args.bin + "/savi_poster_accum > " + itdir + "/" + str(j)
				
				# run cmd and store SGE job id
				myjobid = whichcmd(cmd, args, 1, args.sge, "j_s" + str(args.sampleindex) + "_i" + str(i) + "." +  str(j) + "_post")

				if ( args.sge ):
					print("Your job " + myjobid + " has been submitted")
					sgejobids = myjobid + "," + sgejobids

	if (args.iteration > 1):
		cmd="ln -sf res/" + str(args.iteration) + ' ' + args.outputdir + "/prior"
		whichcmd(cmd, args, 0)
		# if args.verbose: print ("\n# plot prior")
		# plot_prior(args, resdir + "/" + str(args.iteration))

	if ( args.verbose ): print("[END]")

# def plot_prior(args, priorfile):
# 	"""Plot graph of prior"""
# 
# 	cmd = "cat " + priorfile + " | " + args.bin + "/savi_poster_prt"
# 	priorxy = whichcmd(cmd, args, 1)
# 	pylab.xlabel('frequency')
# 	pylab.title("Prior for " + args.name )
# 	pylab.plot(priorxy.split()[0::2], priorxy.split()[1::2])
# 	pylab.savefig(args.outputdir + "/prior.png")

def get_arg():
	"""Get Arguments"""

	prog_description = """Run savi pipeline - PRIOR MAKING STEP.
	You must supply a starting prior file and a compressed (multi-sample) vcf file,
	and tell this script which sample to use (--sampleindex flag) in the vcf 
	as well as which regions (--region flag).
	Notes: This wrapper needs to find the savi binaries in order to function properly. 
	It assumes they are in a directory called savi/ in the same folder where this script resides.
	If this is not the case, you can supply the path.
	You must also supply a starting prior---there are some available to use in the savi/ directory mentioned above---and
	a temporary directory for the script to work in.
	(If you don't, it will make a folder tmp/ in the cwd and work there)
	If you want to run the script in the style of the old savi (not recommended), you can supply a comma-delimited list of (1,q,v,t) files instead of a vcf.
	This script has the following dependencies: tabix and bgzip (available here: http://samtools.sourceforge.net/tabix.shtml);
	vcflib (available here: https://github.com/ekg/vcflib);
	and pylab (http://wiki.scipy.org/PyLab)
	"""

	# http://docs.python.org/2/howto/argparse.html
	parser = argparse.ArgumentParser(description=prog_description)
	parser.add_argument("-v", "--verbose", 						action="store_true",		help="verbose mode. Default: off")
	parser.add_argument("--version", 						action="store_true",		help="print version and exit")
	parser.add_argument("--example", 						action="store_true",		help="print example command and exit")
	# argument as well as pipe-able
	parser.add_argument("-i", "--input",										help="a (multi-sample) vcf file whose sample fields at least have the subfields RBQ ABQ RD AD. Must be bgzipped and indexed with tabix")
	parser.add_argument("-r", "--region",										help="a comma-delimited list of the regions in your vcf you want to build the prior on. If you dont specify this, the prior will run on every line of the vcf file")
	parser.add_argument("-s", "--sampleindex",	type=int,			default=1,			help="index of the sample to use in your multi-sample vcf. Default: 1")
	parser.add_argument("--qvt",										help="If you dont supply a vcf file, you may instead supply a comma-delimited list of 1-qual-vardepth-totdepth files, the way the old scripts used to run. This is not recommended")
	parser.add_argument("-n", "--name",						default="sample",		help="the name of your sample. Default: sample")
	parser.add_argument("-j", "--iteration",	type=int,			default=10,			help="number of iterations. Default: 10")
	parser.add_argument("-p", "--prior",										help="first prior (e.g., \"unif_prior01\")")
	parser.add_argument("-b", "--bin",						default = script_dir + "/savi",	help="directory where the savi binaries reside. Default: [scriptdir]/savi")
	parser.add_argument("-o", "--outputdir",					default = cwdir + "/tmp", 	help="output directory. Default: [cwd]/tmp")
	parser.add_argument("-q", "--sge",						action="store_true",		help="qsub commands with SGE. Default: off")
	parser.add_argument("-l", "--sgelog",						default = cwdir + "/log", 	help="for SGE only. Default: [cwd]/log")
	parser.add_argument("-m", "--sgemem",						default ="4", 			help="for SGE only. sge memory allotment (hrs). Default: 4")
	parser.add_argument("-t", "--sgetime",						default ="4", 			help="for SGE only. sge time allotment (Gigabytes). Default: 4")
	args = parser.parse_args()

	# check for any input
	if ( not len(sys.argv) > 1 ):
		parser.print_help()
		sys.exit(0)

	# print version information 
	if ( args.version ):
		print( os.path.basename(__file__) + ": Version 1.0" )
		sys.exit(0)

	# print example command 
	if ( args.example ):
		my_example = " -v -n test -j 5 -p ../savi/unif_prior01 -o dir_prior -i test.vcf.bgz -r 1:16000-17000,1:18000-18500,1:19000-20000\n\n" + \
			     "This would run 5 iterations of prior building starting with savi/unif_prior01.\n" + \
			     "The raw material would be the first sample in test.vcf.bgz, regions chr1:16000-17000, 18000-18500, and 19000:20000.\n" + \
			     "The output would be saved in dir_prior/"
		print( "# " + os.path.basename(__file__) + my_example )
		sys.exit(0)

	# check for input vcf
	if ( args.input and not os.path.isfile(args.input) ):
		print("can't find input vcf " + args.input)
		sys.exit(1)

	# check for input qvt, build up abspathstr
	if ( args.qvt and not args.qvt == '-' ):
		# store the absolute paths
		abspathstr = []

		for i in args.qvt.split(","):
			if ( os.path.isfile(i) ):
				abspathstr.append(os.path.abspath(i))
			else:
				print("can't find input qvt " + i)
				sys.exit(1)

		# convert arguments to absolute path 
		args.qvt = ",".join(abspathstr)

	if ( not args.input and not args.qvt ):
		print("can't find input")
		sys.exit(1)

	# check for prior
	if ( not os.path.isfile(args.prior) ):
		print("can't find input " + args.prior)
		sys.exit(1)

	# check for binaries
	for j in ("savi_poster", "savi_poster_accum", "savi_poster_merge", "make_qvt"):
		if ( not os.path.isfile(args.bin + "/" + j) ):
			print("can't find " + args.bin + "/" + j)
			sys.exit(1)
	
	if ( args.verbose ):
		print("[BEGIN]")
		print("Arguments: "),
		print(args)

	return args

# -------------------------------------

if __name__ == "__main__":

	main()	
