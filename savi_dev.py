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

	# get arguments dictionary
	args = get_arg()

	# check for errors, check dependencies
	check_error(args)

	# if user requests, run the appropriate step
	for i in args.steps:
		# get the class
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
	parser.add_argument("--region","-r",			default = '',	help="the genomic region to run SAVI on (default: full range) (example: chr1 or chr1:1-50000000)")
	# parser.add_argument("--index",						help="the integer index of the region you want to run SAVI on. By default, 1 refers to the first chromosome in range.txt, 2 to the second, and so on. If you use the partition flag, 1 refers to the first partition, 2 to the second and so on.  (default: 0, which corresponds to the full range)")
	# parser.add_argument("--partition",						help="Number of bases into which to partition the genome.  Use this flag, only if you want to break your genome into regions smaller than the chr length.  If you use this flag, you must also specify an index. For example, \"--partition 50000000 --index 1\" would refer to chr1:1-50000000 for hg19 (default: not used)")
	parser.add_argument("--names",							help="sample names in a comma-delimited list, in the corresponding order of your bam files (default: names are numerical indicies)")
	parser.add_argument("--compsamp","-c",						help="comma-delimited list of colon-delimited indices of samples to compare with savi (default: everything compared to sample 1) (example: 2:1 would compare the second bam file to the first) (example: 2:1,3:1,3:2 would compare the second to the first, the third to the first, and the third to the second)")
	parser.add_argument("--steps",				default="1245",		help="steps to run (default: 1,2,4,5 (i.e., all except prior generation))")
	parser.add_argument("--ann-genome",			default="hg19",		help="name of the SnpEff genome with which to annotate (default: hg19)")
	parser.add_argument("--memory",		type=int,	default=6,		help="the memory for the (SnpEff) Java virtual machine in gigabytes (default: 6)")
	parser.add_argument("--scripts",			default=software,	help="location of scripts dir (directory where this script resides - use this option only if qsub-ing with the Oracle Grid Engine)")
	parser.add_argument("--mindepth",	type=int,	default=10,		help="the min tot read depth required in at least one sample - positions without this wont appear in pileup file (default: 10). Where the filtering occurs: samtools mpileup post-processing")
	parser.add_argument("--minad",		type=int,	default=2,		help="the min alt depth (AD) in at least one sample to output variant (default: 2). Where the filtering occurs: samtools mpileup post-processing")
	parser.add_argument("--mapqual",	type=int,	default=10,		help="skip alignments with mapQ less than this (default: 10). Where the filtering occurs: samtools mpileup")
	parser.add_argument("--maxdepth",	type=int,	default=100000,		help="max per-BAM depth option for samtools mpileup (default: 100000)")
	parser.add_argument("--s1adpp",		type=int,	default=3,		help="for filtered report, require the sample 1 (normal) alt depth per position to be less than this (default: 3) (note: this is NOT sample1 alt depth of the given alt but, rather, at the given position). Where the filtering occurs: generating report.coding.somatic")
	parser.add_argument("--minallelefreq",		type=int,	default=4,	help="Sgt1MAXFREQ (the allele frequency of any sample not including the first one, assumed to be normal) is greater than this (default: 4) Where the filtering occurs: generating the PD.report file.")
	parser.add_argument("--ann-vcf",						help="comma-delimited list of vcfs with which to provide additional annotation (default: none). Where it's used: SnpSift")
	parser.add_argument("--buildprior",			default=software + "/bin/prior_unif01",		help="starting input prior when building the prior if step 3 is invoked (default: bin/prior_unif01)")
	parser.add_argument("--prior",				default=software + "/bin/prior_diploid01",	help="prior to use if step 3 is not run (default: bin/prior_diploid01)")
	parser.add_argument("--prioriterations",		default="10",		help="the number of iterations for the prior build, if step 3 is run (default: 10)")
	parser.add_argument("--presence",			default="1e-6",		help="the SAVI presence posterior (default: 1e-6). Where it's used: step 4")
	parser.add_argument("--conf",				default="1e-5",		help="the SAVI conf (default: 1e-5). Where it's used: step 4")
	parser.add_argument("--precision",	type=int,	default=0,		help="the SAVI precision (default: 0). Where it's used: step 4")
	parser.add_argument("--rnabams",						help="comma-delimited list of rna bam files (.bai indices should be present)")
	# parser.add_argument("--nofilter",	action="store_true",			help="do not use SAVI comparison filter (default: off) (you should use this option if NOT doing comparisons - e.g., if you only have one sample)")
	parser.add_argument("--nopv4",		action="store_true",			help="do not run bcftools to compute PV4 (default: off). Where it's used: step 5")
	parser.add_argument("--noindeldepth",	action="store_true",			help="do not include include indel reads in total depth count (SDP) (default: off) (note: asteriks in mpileup are included irrespective of this flag). Where it's used: step 2")
	parser.add_argument("--rdplusad",	action="store_true",			help="use reference-agreeing reads plus alternate-calling reads (RD+AD) rather than total depth (SDP) as input to savi (default: off). Where it's used: step 2")
	# parser.add_argument("--hybrid",	action="store_true",			help="as input to savi, use reference-agreeing reads plus alternate-calling reads (RD+AD) for first sample (normal) and SDP for other samples (default: off) (note: this flag changes read depths on positions where there are multiallelic variants)")
	parser.add_argument("--index",				default="0",		help="an index used in naming of output files (default: 0)")
	parser.add_argument("--noncoding",	action="store_true",			help="use snpEff to find all transcripts, not just only protein transcripts (default: off). Where it's used: step 5")
	parser.add_argument("--noclean",	action="store_true",			help="do not delete temporary intermediate files (default: off)")
	parser.add_argument("--debug",		action="store_true",			help="echo commands (default: off)")

	args = parser.parse_args()

	# if no args, print help
	if ( not len(sys.argv) > 1 ):
		parser.print_help()
		sys.exit(0)

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
	# generate default prior string
	vars(args)['priorstring'] = generate_priorstr(args.compsamp, args.prior)

	# need to fix paths if users passed software path manually

	# print args
	print(args)
	print

	return args

# -------------------------------------

def check_error(args):
	"""Check for errors, check dependencies """

	# check step string only goes from 1 to 5 and is all digits 
	if not args.steps.isdigit():
		print("[ERROR] invalid step string")
		sys.exit(1)
	if set(list(args.steps)).intersection(set(list('67890'))):
		print("[ERROR] invalid step string (only steps 1 thro 5 exist)")
		sys.exit(1)

	# required programs
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
	return ",".join(["%s:%s" % (k, priorpath) for k in get_uniq_samples(compsamp)])

# -------------------------------------

def run_cmd(cmd, bool_verbose, bool_getstdout):
	"""Run system command"""

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
	"""check for existence and nonzero-ness of output file"""

	if (os.path.isfile(myfile)):
		if (os.path.getsize(myfile) == 0):
			print(myfile + " is empty. Exiting")
			sys.exit(1)
	else:
		print("Can't find " + myfile + ". Exiting.")
		sys.exit(1)

# -------------------------------------

class Step(object):
	"""A parent step class from which step children inherit"""

	def __init__ (self, args, step_index):
		"""initialize step object"""
		# set arguments dictionary
		self.args = args 
		# set step index
		self.step_index = step_index 

	def run(self):
		"""the run method, meant to be overridden by the children"""
		print('[STEP ' + self.step_index + ']')

# -------------------------------------

class Step1(Step):
	"""A step1 object whose run method converts bams to mpileup"""

	def checkok(self):
		"""check if everything's ok - i.e., inputs exist"""

		if not self.args.bams:
			print("[ERROR] No bam files given as input")
			sys.exit(1)

		if not self.args.ref:
			print("[ERROR] No reference given as input")
			sys.exit(1)

	# override parent's run method
	def run(self):
		"""The run method calls shell(system) commands to do step1 - namely, it calls samtools mpileup to convert bam to mpileup"""

		# define input and output attributes
		self.input = self.args.bams
		self.output = self.args.outputdir + "/" + "tmp_mpile." + self.args.index + ".txt"

		# run parent's run method
		super(Step1, self).run()
		# check if there are errors
		self.checkok()
		# define pileup file
		# note the -d flag: "max per-BAM depth to avoid excessive memory usage [250]." 
		# The default of 250 is way too low, so we jack it up
		pileupflag="-A -B -q {} -d {} -L {} -f {}".format(self.args.mapqual, self.args.maxdepth, self.args.maxdepth, self.args.ref)

		# define command to run

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
		run_cmd(mycmd, self.args.debug, 1)

		# check if output file nonzero size
		check_file_exists_and_nonzero(self.output)

# -------------------------------------

class Step2(Step):
	"""A step2 object whose run method converts mpileup to vcf"""

	# override parent's run method
	def run(self):
		"""The run method calls shell(system) commands to do step2 - namely, it converts mpileup to vcf"""

		# define input and output attributes
		self.input = self.args.outputdir + "/" + "tmp_mpile." + self.args.index + ".txt"
		self.output = self.args.outputdir + "/" + self.args.index + ".vcf.bgz"

		# run parent's run method
		super(Step2, self).run()

		# check if input file nonzero size
		check_file_exists_and_nonzero(self.input)

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
			run_cmd(mycmd, self.args.debug, 1)

			# zip and tabix all variants file, then remove unzipped file
			mycmd = "cat {} | bgzip > {}.bgz; rm {}; tabix -p vcf {}.bgz".format(allvariants, allvariants, allvariants, allvariants)
			run_cmd(mycmd, self.args.debug, 1)

		# otherwise, just generate variants file
		else:
			mycmd = "cat {} | {}/pileup2multiallele_vcf {} | bgzip > {}".format(self.input, self.args.bin, oscanflag, self.output)
			run_cmd(mycmd, self.args.debug, 1)

		# tabix variants
		mycmd = "tabix -p vcf {}".format(self.output)
		run_cmd(mycmd, self.args.debug, 1)

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

	def run(self):
		"""The run method calls shell(system) commands to do step3 - namely, compute the prior"""

		# run parent's run method
		super(Step3, self).run()

		# define input and output attributes
		self.input = self.args.outputdir + "/" + self.args.index + ".all.vcf.bgz"
		# self.output = 

		# check if input file nonzero size
		check_file_exists_and_nonzero(self.input)

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
			myout = run_cmd(mycmd, self.args.debug, 1)
			if (self.args.debug): print(myout)

		# if run step3, create prior string for priors specified in compsamp
		# if compsamp="2:1", it should look like this - 1:${outputdir}/savi/prior1/prior,2:${outputdir}/savi/prior2/prior
		# priorstring=$( for i in $( echo ${compsamp} | sed 's|,|\n|g; s|:|\n|g' | sort -u ); do echo -n "${i}:${outputdir}/savi/prior${i}/prior,"; done | sed 's|,$||'; echo )
		# if not step3, use diploid prior by default
		# priorstring=$( for i in $( echo ${compsamp} | sed 's|,|\n|g; s|:|\n|g' | sort -u ); do echo -n "${i}:${inputprior},"; done | sed 's|,$||'; echo )

# -------------------------------------

class Step4(Step):
	"""A step4 object whose run method runs savi proper"""

	# override parent's run method
	def run(self):
		"""The run method calls shell(system) commands to do step4 - namely, call savi proper"""

		# define input and output attributes
		self.input = self.args.outputdir + "/" + self.args.index + ".vcf.bgz"
		# self.output = 

		# run parent's run method
		super(Step4, self).run()

		# check if input file nonzero size
		check_file_exists_and_nonzero(self.input)

		# define a directory for reports
		reportdir = self.args.outputdir + "/report"

		# if dir doesn't exist, create it
		if not os.path.exists(reportdir): os.makedirs(reportdir)

		# add this key-value pair to the args dict
		vars(self.args)['report'] = reportdir

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
		    " --outputdir " + reportdir

		# run it
		myout = run_cmd(mycmd, self.args.debug, 1)
		if (self.args.debug): print(myout)

		# run_savi.py output files:
		# freqsavi.vcf.bgz - adds presence frequencies to all present variants 
		# finalsavi.vcf.bgz - add frequency deltas for sample comparisons to all present variants 
		# finalfilter.vcf - filter all present mutations for somatic variants

# -------------------------------------

class Step5(Step):
	"""A step5 object whose run method annotates with SnpEff and generates report"""

	def run(self):
		"""The run method calls shell(system) commands to do step5 - namely, annotate the vcf"""
		pass

# -------------------------------------

if __name__ == "__main__":

	main()	
