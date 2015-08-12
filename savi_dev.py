#!/usr/bin/env python

"""
	SAVI
	~~~~~~
	SAVI is a program for calling variants in high-throughput sequencing data, 
	particularly paired tumor-normal samples, in 5 distinct steps.
"""

__author__ = "Oliver"
__version__ = "Revision: 2.0"
__date__ = "Date: 07-2015"

import argparse
import sys 
import re
import subprocess
import os
from distutils import spawn

# -------------------------------------

def main():
	"""Main block"""

	# print(globals())

	# get arguments dictionary
	(args, parser) = get_arg()

	# check for errors, check dependencies
	check_error(args, parser)

	# if user requests, run the appropriate step
	for i in args.steps:
		# get the class from the classes in the current script
		myclass = globals()['Step' + i]
		# instantiate a step obj for the class
		mystep = myclass(args, i)
		# execute the run method
		mystep.run()

# -------------------------------------

def get_arg():
	"""Get Arguments"""

	prog_description = """SAVI is a program for calling variants in high-throughput sequencing data, particularly paired tumor-normal samples. In bioinformatics, \"calling variants\" can mean two different things: (1) simply enumerating all differences between some sequence data and a reference; and (2) determining which of those differences are significant and not likely to be error. SAVI does the later, while a program like Samtools mpileup will do the former. In practice, SAVI is a way of sifting through large amounts of data to pull out the significant mutations using a Bayesian probability model. A common use case is identifying deleterious mutations in cancer, given normal and tumor sequence data---an often-encountered problem in bioinformatics. The output of SAVI is a list of candidate genomic alterations each annotated with a probability of how likely it is to be real. SAVI works with standard bioinformatic file formats. For a complete description of how to use this software, including dependencies and usage examples, see https://github.com/RabadanLab/SAVI"""

	# directory where this script resides 			
	software = os.path.dirname(os.path.realpath(__file__))
	# cwd
	cwdir = os.getcwd()

	# http://docs.python.org/2/howto/argparse.html
	parser = argparse.ArgumentParser(description=prog_description)

	parser.add_argument("--bams","-b",						help="comma-delimited list of bam files (by convention list the normal sample first, as in: normal.bam,tumor.bam) (.bai indices should be present)")
	parser.add_argument("--ref",							help="reference fasta file with faidx index in the same directory")
	parser.add_argument("--outputdir","-o",			default = cwdir,	help="the output directory (default: cwd)")
	parser.add_argument("--region","-r",			default = '',		help="the genomic region to run SAVI on (default: full range) (example: chr1 or chr1:1-50000000)")
	parser.add_argument("--names",							help="sample names in a comma-delimited list, in the corresponding order of your bam files (default: names are numerical indicies)")
	parser.add_argument("--compsamp","-c",						help="comma-delimited list of colon-delimited indices of samples to compare with savi (default: everything compared to sample 1) (example: 2:1 would compare the second bam file to the first) (example: 2:1,3:1,3:2 would compare the second to the first, the third to the first, and the third to the second)")
	parser.add_argument("--steps",				default="1245",		help="steps to run (default: 1,2,4,5 (i.e., all except prior generation))")
	parser.add_argument("--ann",				default="hg19",		help="name of the SnpEff genome with which to annotate (default: hg19)")
	parser.add_argument("--memory",				default="4",		help="the memory for the (SnpEff) Java virtual machine in gigabytes (default: 6)")
	parser.add_argument("--scripts",			default=software,	help="location of scripts dir (directory where this script resides - use this option only if qsub-ing with the Oracle Grid Engine)")
	parser.add_argument("--mindepth",	type=int,	default=10,		help="the min tot read depth required in at least one sample - positions without this wont appear in pileup file (default: 10). Where the filtering occurs: samtools mpileup post-processing")
	parser.add_argument("--minad",		type=int,	default=2,		help="the min alt depth (AD) in at least one sample to output variant (default: 2). Where the filtering occurs: samtools mpileup post-processing")
	parser.add_argument("--mapqual",	type=int,	default=10,		help="skip alignments with mapQ less than this (default: 10). Where the filtering occurs: samtools mpileup")
	parser.add_argument("--maxdepth",	type=int,	default=100000,		help="max per-BAM depth option for samtools mpileup (default: 100000)")
	parser.add_argument("--s1adpp",		type=int,	default=3,		help="for filtered report, require the sample 1 (normal) alt depth per position to be less than this (default: 3) (note: this is NOT sample1 alt depth of the given alt but, rather, at the given position). Where the filtering occurs: generating report.coding.somatic")
	parser.add_argument("--minallelefreq",		type=int,	default=4,	help="Sgt1MAXFREQ (the allele frequency of any sample not including the first one, assumed to be normal) is greater than this (default: 4) Where the filtering occurs: generating the PD.report file.")
	parser.add_argument("--annvcf",							help="comma-delimited list of vcfs with which to provide additional annotation (default: none). Where it's used: SnpSift")
	parser.add_argument("--buildprior",			default=software + "/bin/prior_unif01",		help="starting input prior when building the prior if step 3 is invoked (default: bin/prior_unif01)")
	parser.add_argument("--prior",				default=software + "/bin/prior_diploid01",	help="prior to use if step 3 is not run (default: bin/prior_diploid01)")
	parser.add_argument("--prioriterations",		default="10",		help="the number of iterations for the prior build, if step 3 is run (default: 10)")
	parser.add_argument("--presence",			default="1e-6",		help="the SAVI presence posterior (default: 1e-6). Where it's used: step 4")
	parser.add_argument("--conf",				default="1e-5",		help="the SAVI conf (default: 1e-5). Where it's used: step 4")
	parser.add_argument("--precision",	type=int,	default=0,		help="the SAVI precision (default: 0). Where it's used: step 4")
	parser.add_argument("--rnabams",						help="comma-delimited list of rna bam files (.bai indices should be present)")
	parser.add_argument("--nopv4",		action="store_true",			help="do not run bcftools to compute PV4 (default: off). Where it's used: step 5")
	parser.add_argument("--noindeldepth",	action="store_true",			help="do not include include indel reads in total depth count (SDP) (default: off) (note: asteriks in mpileup are included irrespective of this flag). Where it's used: step 2")
	parser.add_argument("--rdplusad",	action="store_true",			help="use reference-agreeing reads plus alternate-calling reads (RD+AD) rather than total depth (SDP) as input to savi (default: off). Where it's used: step 2")
	# parser.add_argument("--hybrid",	action="store_true",			help="as input to savi, use reference-agreeing reads plus alternate-calling reads (RD+AD) for first sample (normal) and SDP for other samples (default: off) (note: this flag changes read depths on positions where there are multiallelic variants)")
	parser.add_argument("--index",				default="0",		help="an index used in naming of output files (default: 0)")
	parser.add_argument("--noncoding",	action="store_true",			help="use snpEff to find all transcripts, not just only protein transcripts (default: off). Where it's used: step 5")
	parser.add_argument("--noclean",	action="store_true",			help="do not delete temporary intermediate files (default: off)")
	parser.add_argument("--verbose","-v",	action="store_true",			help="echo commands (default: off)")
	parser.add_argument("--superverbose",	action="store_true",			help="echo commands (default: off)")

	args = parser.parse_args()

	# automatically modify the argument dict:

	# get the number of samples
	numsamp = 0
	if args.bams:
		numsamp = len(args.bams.split(",")) 

	# if user passes scripts dir, change paths of priors (default and build)
	if args.scripts:
		# change default if user didn't supply his own input
		if args.buildprior == software + "/bin/prior_unif01":
			args.buildprior = args.scripts + "/bin/prior_unif01"
		if args.prior == software + "/bin/prior_diploid01":
			args.prior = args.scripts + "/bin/prior_diploid01"

	# add key-value pairs to the args dict
	vars(args)['numsamp'] = numsamp 
	vars(args)['cwd'] = cwdir
	vars(args)['bin'] = args.scripts + "/bin" 
	# if user didn't define compsamp, set it
	if not args.compsamp: vars(args)['compsamp'] = generate_compsamp(numsamp)
	# generate default prior string if user not running step3
	if not "3" in args.steps:
		vars(args)['priorstring'] = generate_priorstr(args.compsamp, args.prior)
	# define a directory for reports
	vars(args)['reportdir'] = args.outputdir + "/report"

	# print args
	print(args)
	print

	return args, parser

# -------------------------------------

def check_error(args, parser):
	"""Check for errors, check dependencies """

	# if no args, print help
	if ( not len(sys.argv) > 1 ):
		parser.print_help()
		sys.exit(0)

	# unpaired?
	if (args.numsamp == 1 ):
		print("[ALERT] running savi with an unpaired sample - this is not optimal")

	# check step string only goes from 1 to 5 and is all digits 
	if not args.steps.isdigit():
		print("[ERROR] invalid step string")
		sys.exit(1)
	if set(list(args.steps)).intersection(set(list('67890'))):
		print("[ERROR] invalid step string (only steps 1 thro 5 exist)")
		sys.exit(1)

	# define required programs
	required_programs = ("java", "python", "samtools", "snpEff.jar", "bgzip", "tabix", "vcffilter", "bcftools")
	required_savi_programs = ("pileup2multiallele_vcf", "add_to_info", "make_qvt", "savi_poster", "savi_conf", "savi_comp", "savi_poster_accum", "savi_poster_merge", "savi_poster_prt", "savi_unif_prior", "savi_txt2prior")

	# check for required programs 
	for j in required_savi_programs:
		check_path(args.bin + "/" + j)

	for j in required_programs:
		if (not spawn.find_executable(j)):
			print("[ERROR] Can't find " + j + ". Please add it to your PATH")
			sys.exit(1)

	# check for existence of bam files
	if args.bams:
		for j in args.bams.split(","):
			check_path(j)
			check_path(j + '.bai')

	# check for existence of reference
	if args.ref:
		check_path(args.ref)
		check_path(args.ref + '.fai')

	# check for existence of annotating vcf
	if args.annvcf:
		for i in args.annvcf.split(","):
			check_path(i)

	# if step1 called, but no bam files supplied, throw error
	if '1' in args.steps and not args.bams:
		print("[ERROR] No bam files given as input")
		sys.exit(1)

	# need to check versions
	# need to print times

# -------------------------------------

def check_path(myfile):
	"""Check for the existence of a file"""

	# expanduser handle tilde
	if (not os.path.isfile(os.path.expanduser(myfile))):
		print("[ERROR] Can't find the file " + myfile)
		sys.exit(1)

# -------------------------------------

def remove_files(myfile, boolverbose):
	"""Remove files"""

	# command: samtools plus awk
	mycmd = "rm -f " + myfile
	# run it
	run_cmd(mycmd, boolverbose, 1)

# -------------------------------------

def generate_compsamp(numsamp):
	"""Generate sample comparison string, which determines which samples are compared to which"""

	# by default, everything is compared to sample 1, which is assumed to be 'normal' (i.e., not tumor)
	if (numsamp == 1):
		return '1'
	elif (numsamp > 1):
		return ",".join(["%s:1" % k for k in range(2,numsamp+1)]) 
	else:
		return ''

# -------------------------------------

def get_uniq_samples(compsamp):
	"""Return the list of uniq samples from the sample comparison string"""

	return sorted(set(compsamp.replace(':', ',').split(',')))

# -------------------------------------

def generate_priorstr(compsamp, priorpath):
	"""Generate prior string for run_savi.py"""

	# if run step3, create prior string for priors specified in compsamp
	# if compsamp="2:1", it should look like this - 1:${outputdir}/savi/prior1/prior,2:${outputdir}/savi/prior2/prior
	# if not step3, use diploid prior by default
	if not compsamp:
		return "1:" + priorpath
	else:
		return ",".join(["%s:%s" % (k, priorpath) for k in get_uniq_samples(compsamp)])

# -------------------------------------

def run_cmd(cmd, bool_verbose, bool_getstdout):
	"""Run a system (i.e., shell) command"""

	# if verbose, print command
	if (bool_verbose): 
		print("[command] " + cmd)

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	proc.wait() 
	# print return code
	# print(proc.returncode) 
	# print stdout stderr tuple
	# proc.communicate()

	(stdout, stderr) =  proc.communicate()

	# if error, print it
	if stderr:
		print("[stderror] " + stderr),

	# return stdout
	if (bool_getstdout): 
		return stdout.rstrip()
	else:
		return "0" # note: this must return a str

# -------------------------------------

def check_file_exists_and_nonzero(myfile):
	"""Check for the existence and nonzero-ness of a file"""

	# loop through comma-delimited list of files
	for i in myfile.split(","):
		# if (os.path.isfile(i)):
		if (os.path.isfile(os.path.expanduser(i))):
			if (os.path.getsize(i) == 0):
				print(i + " is empty. Exiting")
				sys.exit(1)
		else:
			print("Can't find " + i + ". Exiting.")
			sys.exit(1)

# -------------------------------------

class Step(object):
	"""A parent step class from which step children inherit"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""
		# set arguments dictionary
		self.args = args 
		# set step index
		self.step_index = step_index 

		# list of intermediate files (to be cleaned)
		self.intermediates=[]

	def set_input(self, myfile):
		"""Set the input file for the step"""
		self.input = myfile

		# check if input file nonzero size
		check_file_exists_and_nonzero(self.input)

	def set_output(self, myfile):
		"""Set the output file for the step"""
		self.output = myfile

	def set_descrip(self, mydescrip):
		"""Set the description"""
		self.description = mydescrip 

	def run(self):
		"""The run method, meant to be overridden by the children"""
		print('[STEP ' + self.step_index + '] ' + self.description + '\n')

	def cleanup(self):
		"""Clean up tmp and intermediate files"""
		for i in self.intermediates:
			remove_files(i, self.args.verbose)

# -------------------------------------

class Step1(Step):
	"""A step1 object whose run method converts bams to mpileup"""

	def __init__ (self, args, step_index):
		"""Step1 constructor"""

		# parent's constructor
		Step.__init__(self, args, step_index)

		# define description as well as input and output attributes
		self.set_descrip("Convert bams to mpileup and, if not building custom prior, filter for variants only")
		self.set_input(self.args.bams)
		self.set_output(self.args.outputdir + "/" + "tmp_mpile." + self.args.index + ".txt")

	# override parent's run method
	def run(self):
		"""The run method calls shell(system) commands to do step1 - namely, it calls samtools mpileup to convert bam to mpileup"""

		# run parent's run method
		super(Step1, self).run()

		# define pileup file
		# note the -d flag: "max per-BAM depth to avoid excessive memory usage [250]." 
		# The default of 250 is way too low, so we jack it up
		pileupflag="-A -B -q {} -d {} -L {} -f {}".format(self.args.mapqual, self.args.maxdepth, self.args.maxdepth, self.args.ref)

		# define command to run

		# **TO DO**: port this awk to Python

		# first a line of awk: we want only lines where normal has depth > cutoff
		# and at least one of the tumor samples has depth > cutoff
		# and we want only variants
		awkcmd = "awk -v num=" + str(self.args.numsamp) + " -v cutoff=" + str(self.args.mindepth) + " -v minad=" + str(self.args.minad) + """ '{
			# reduce pileup by only printing variants

			myflag=0; # flag to print line

			# loop thro samples
			for (i = 1; i <= num-1; i++) 
			{
				# only consider if sample depth greater than or eq to cutoff
				if ($(4+i*3) >= cutoff) 
				{
					# get read string
					mystr=$(5+i*3); 
					# eliminate read starts with qual scores ACTG (e.g., stuff like ^A ^C etc)
					gsub(/\^[ACTGNactgn+-]/,"",mystr); 
					# count ACTGN mismatches
					altdepth = gsub(/[ACTGNactgn]/,"",mystr); 
					# if sufficient number of mismatches, flag line to print
					if (altdepth >= minad) 
					{
						myflag=1; break;
					}
				}
			}

			# if only one sample, print if depth >= cutoff
			if (num==1 && $4 >= cutoff) 
			{
				print;
			}

			# if multiple samples, print if sample one has depth AND myflag
			else if (num>1 && $4 >= cutoff && myflag) 
			{
				print;
			}
		}'"""

		# **TO DO**: port this awk to Python

		# if make prior step, modify awk command to generate all positions in the mpileup, variants and non-variants alike
		if "3" in self.args.steps:
			awkcmd = "awk -v num=" + str(self.args.numsamp) + " -v cutoff=" + str(self.args.mindepth) + """ '{
				true=0; 
				for (i = 1; i <= num-1; i++) 
				{
					if ($(4+i*3) >= cutoff) {true=1}
				}; 
				if (num==1 && $4 >= cutoff) 
				{
					print;
				}
				else if (num>1 && $4 >= cutoff && true)
				{
					print;
				}
			}'"""

		# command: samtools plus awk
		mycmd = "samtools mpileup {} {} {} | {} > {}".format(pileupflag, self.args.region, self.args.bams.replace(',', ' '), awkcmd, self.output)
		# run it
		run_cmd(mycmd, self.args.verbose, 1)

		# check if output file nonzero size
		check_file_exists_and_nonzero(self.output)

# -------------------------------------

class Step2(Step):
	"""A step2 object whose run method converts mpileup to vcf"""

	def __init__ (self, args, step_index):
		"""Step2 constructor"""

		# parent's constructor
		Step.__init__(self, args, step_index)

		# define description as well as input and output attributes
		self.set_descrip("Convert mpileup to vcf")
		self.set_input(self.args.outputdir + "/" + "tmp_mpile." + self.args.index + ".txt")
		self.set_output(self.args.outputdir + "/" + self.args.index + ".vcf.bgz")

	# override parent's run method
	def run(self):
		"""The run method calls shell(system) commands to do step2 - namely, it converts mpileup to vcf"""

		# run parent's run method
		super(Step2, self).run()

		# define command to run

		# first, define flag for pileup2multiallele_vcf (nicknamed oscan)
		# addresses a problem that samtools mpileup does not include indels in the total read count (with the exception of *s)
		# so add indel depths to SDP
		oscanflag="--cutoff " + str(self.args.minad) + " --header --includeindels"

		# if requested, don't include indel reads in total depth count
		if (self.args.noindeldepth):
			oscanflag="--cutoff " + str(self.args.minad) + " --header"

		# if make prior, need to generate all position vcf (variants and non-variants)
		if "3" in self.args.steps:

			# add --all to flag
			oscanflag += " --all"

			# command: pileup2multiallele_vcf
			awkcmd = """awk '{ if ($0 ~ /^#/ || toupper($4) != toupper($5)) {print}}'"""
			allvariants = self.args.outputdir + "/" + self.args.index + ".all.vcf"
			mycmd = "cat {} | {}/pileup2multiallele_vcf {} | tee {} | {} | bgzip > {}".format(self.input, self.args.bin, oscanflag, allvariants, awkcmd, self.output)
			# run it
			run_cmd(mycmd, self.args.verbose, 1)

			# zip and tabix all variants file, then remove unzipped file
			mycmd = "cat {} | bgzip > {}.bgz; rm {}; tabix -p vcf {}.bgz".format(allvariants, allvariants, allvariants, allvariants)
			run_cmd(mycmd, self.args.verbose, 1)

		# otherwise, just generate variants file
		else:
			mycmd = "cat {} | {}/pileup2multiallele_vcf {} | bgzip > {}".format(self.input, self.args.bin, oscanflag, self.output)
			run_cmd(mycmd, self.args.verbose, 1)

		# tabix variants
		mycmd = "tabix -p vcf {}".format(self.output)
		run_cmd(mycmd, self.args.verbose, 1)

		# if empty
		# if [ $( zcat ${outputdir}/${SGE_TASK_ID}.vcf.bgz | sed '/^#/d' | head | wc -l ) == 0 ]; then
		# 	echo "[HALT SCRIPT] vcf file is empty"
		# 	exit 0;
		# fi

		# if [ $cleanup -eq 1 ]; then 
		# 	rm -f ${outputdir}/tmp_mpile.${SGE_TASK_ID}.*
		# fi

# -------------------------------------

class Step3(Step):
	"""A step3 object whose run method computes prior"""

	def __init__ (self, args, step_index):
		"""Step3 constructor"""

		# parent's constructor
		Step.__init__(self, args, step_index)

		# define description as well as input and output attributes
		self.set_descrip("Construct the prior for each sample")
		self.set_input(self.args.outputdir + "/" + self.args.index + ".all.vcf.bgz")

	def run(self):
		"""The run method calls shell(system) commands to do step3 - namely, compute the prior"""

		# run parent's run method
		super(Step3, self).run()

		# make prior string
		vars(self.args)['priorstring'] = "1:" + self.args.outputdir + "/priors/prior1/prior"

		# loop through unique samples
		for i in get_uniq_samples(self.args.compsamp):

			# define a directory for prior
			priordir = self.args.outputdir + "/priors/prior" + i

			# if dir doesn't exist, create it
			if not os.path.exists(priordir): os.makedirs(priordir)

			# command: make_prior
			mycmd = self.args.scripts + "/make_prior.py --verbose --name s" + i + \
			    " --input " + self.input + \
			    " --iteration " + self.args.prioriterations + \
			    " --sampleindex " + i + \
			    " --prior " + self.args.buildprior + \
			    " --outputdir " + priordir 

			# run it
			myout = run_cmd(mycmd, self.args.verbose, 1)
			if (self.args.superverbose): print(myout)

			# update prior str
			if not i == "1":
				vars(self.args)['priorstring'] += "," + i + ":" + self.args.outputdir + "/priors/prior" + i + "/prior"

		if self.args.verbose: print("[priorstring] " + vars(self.args)['priorstring']) 

# -------------------------------------

class Step4(Step):
	"""A step4 object whose run method runs savi proper"""

	def __init__ (self, args, step_index):
		"""Step4 constructor"""

		# parent's constructor
		Step.__init__(self, args, step_index)

		# define description as well as input and output attributes
		self.set_descrip("Add SAVI statistics to the vcf")
		self.set_input(self.args.outputdir + "/" + self.args.index + ".vcf.bgz")

	# override parent's run method
	def run(self):
		"""The run method calls shell(system) commands to do step4 - namely, call savi proper"""

		# run parent's run method
		super(Step4, self).run()

		# if report dir doesn't exist, create it
		if not os.path.exists(self.args.reportdir): os.makedirs(self.args.reportdir)

		# define command to run

		# first define flag
		saviflag=""
		# if user requests, use AD+RD as tot depth
		if (self.args.rdplusad):
			saviflag="--rdplusad"

		# keep freq file is off by default but turns on if no sample comparisions
		keepfreqfile = "0"
		if not ":" in self.args.compsamp: keepfreqfile = "1"

		# command: run_savi 
		mycmd = self.args.scripts + "/run_savi.py --verbose " + saviflag + \
		    " --input " + self.input + \
		    " --name savi_" + self.args.index + \
		    " --sample " + self.args.compsamp + \
		    " --prior " + self.args.priorstring + \
		    " --saviprecision " + str(self.args.precision) + \
		    " --savipresent " + self.args.presence + \
		    " --saviconf " + self.args.conf + \
		    " --keepfreqfile " + keepfreqfile + \
		    " --outputdir " + self.args.reportdir 

		# run it
		myout = run_cmd(mycmd, self.args.verbose, 1)
		if (self.args.superverbose): print(myout)
		
		# run_savi.py output files:
		# freqsavi.vcf.bgz - adds presence frequencies to all present variants 
		# finalsavi.vcf.bgz - add frequency deltas for sample comparisons to all present variants 
		# finalfilter.vcf - filter all present mutations for somatic variants

		# need to unzip (a bit awkward)
		# if no savi comparisons (e.g., unpaired sample)
		if (int(keepfreqfile)):
			mycmd = "gunzip -S .bgz " + self.args.reportdir + "/freqsavi.vcf.bgz"
		else:
			mycmd = "gunzip -S .bgz " + self.args.reportdir + "/finalsavi.vcf.bgz"

		myout = run_cmd(mycmd, self.args.verbose, 1)

# -------------------------------------

class Step5(Step):
	"""A step5 object whose run method annotates with SnpEff and generates report"""

	def __init__ (self, args, step_index):
		"""Step5 constructor"""

		# parent's constructor
		Step.__init__(self, args, step_index)

		# define description as well as input and output attributes
		self.set_descrip("Add annotation to the vcf")
		self.set_input(self.args.reportdir + "/finalsavi.vcf")

		# set output prefix, output reports
		self.outprefix = self.args.reportdir + "/tmp_" + self.args.index

		# full vcf report 
		self.vcf_all = self.args.reportdir + "/report.all.vcf"
		# vcf report, filtered for coding regions 
		self.vcf_coding = self.args.reportdir + "/report.coding.vcf"
		# vcf report, filtered for coding regions + somatic variants
		self.vcf_somatic = self.args.reportdir + "/report.coding.somatic.vcf"

		# full tsv report
		self.report_all = self.args.reportdir + "/report.unfiltered.txt"

	def runEff(self):
		"""Run SnpEff and SnpSift"""

		# TO DO: implement function to get rid of Ns before annotation because SnpEff messes this up 
		# (see bug report: https://github.com/pcingola/SnpEff/issues/54)

		# run parent's run method
		super(Step5, self).run()

		# define SnpEff command to run
		mycmd = "" 

		effopts = " -noLog -noHgvs -q -formatEff -noStats -lof -canon"

		# if not noncoding option, annotate only somatic variants (not all present)
		if not self.args.noncoding: 
			effopts += " -onlyProtein -no-downstream -no-intergenic -no-upstream -no-utr"

		# get paths
		effpath = spawn.find_executable('snpEff.jar')
		siftpath = spawn.find_executable('SnpSift.jar')
		configpath = spawn.find_executable('snpEff.config')

		# set output file
		self.output = self.outprefix + ".eff_0.all.vcf"

		mycmd = "java -Xmx" + self.args.memory + "G" + \
			" -jar " + effpath + \
			" ann " + self.args.ann + \
			effopts + \
			" -c " + configpath + \
			" " + self.input + " > " + self.output

		myout = run_cmd(mycmd, self.args.verbose, 1)
		if (self.args.superverbose): print(myout)

		# define SnpSift commands to run

		# loop through all annotating vcfs
		if self.args.annvcf:
			for j,i in enumerate(self.args.annvcf.split(",")):
				# set output file
				self.output = self.outprefix + ".eff_" + str(j+1) + ".all.vcf"

				mycmd = "java -Xmx" + self.args.memory + "G" + \
					" -jar " + siftpath + \
					" annotate " + i + \
					" -noLog " + \
					self.outprefix + ".eff_" + str(j) + ".all.vcf > " + self.output

				myout = run_cmd(mycmd, self.args.verbose, 1)
				if (self.args.superverbose): print(myout)

		# rename final iteration of vcf to report
		os.rename(self.output, self.vcf_all)

	def filterCoding(self, report_in, report_out):
		"""Filter for the coding region - 
		get useful variants and variants in the coding region. The idea is to capture these SnpEff features:
		inframe_deletion *inframe_insertion frameshift_variant* initiator_codon_variant missense_variant* synonymous_variant*
		splice_acceptor* splice_donor* splice_region* start_lost* start_gained* stop_gained* stop_lost* stop_retained_variant*"""

		# filter for coding region
		# desired features
		featurelist = ["inframe", "frameshift", "synonymous_variant", "missense_variant", "splice_acceptor", "splice_donor", "splice_region", "start", "stop"]
		with open(report_out, "w") as g:
			with open(report_in, 'r') as f:
				for line in f:
					if (line.startswith("#") or any(k in line.split()[7] for k in featurelist)):
						g.write(line)

	def filterCodingSomatic(self, report_in, report_out):
		"""Filter for the coding + somatic variants"""

		# filtering: 
		# S1 AD PP < cutoff so discard variants with normal alt depth above a certain threshhold
		# discard variants with strand bias
		# discard variants in meganormals
		# discard variants in COMMON dbSNP
		# take out strand bias filter
		# vcffilter -f "S1ADPP < ${s1adpp}" ${outputdir2}/tmp_${SGE_TASK_ID}.eff_${myindex}.coding.vcf | awk '$0 ~ /^#/ || ( $0 !~ /MEGANORMAL_ID/ && $0 !~ /Meganormal186GBM/ && $0 !~ /COMMON\=1/)' > ${outputdir2}/tmp_${SGE_TASK_ID}.eff_${myindex}.coding.somatic.vcf

		pass

	def filterCodingPD(self, report_in, report_out):
		"""Filter for the PD lower bound greater than zero"""

		pass

	def runPostEffProcessing(self):
		"""Do the Post-SnpEff Processing - Filter vcf and convert to human readable report"""

		# convert vcf to tsv (tab-delimited report)
		mycmd = "cat " + self.vcf_all + " | " + self.args.bin + "/vcf2fullreport.py > " + self.report_all
		run_cmd(mycmd, self.args.verbose, 1)

		# filter for coding region
		self.filterCoding(self.vcf_all, self.vcf_coding)

	def run(self):
		"""The run method calls shell(system) commands to do step5 - namely, annotate the vcf"""

		# run SnpEff
		self.runEff()
		# run post-SnpEff filtering, processing
		self.runPostEffProcessing()

# -------------------------------------

if __name__ == "__main__":

	main()
