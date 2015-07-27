#!/usr/bin/env python

"""
	SAVI
	~~~~~~
	SAVI is a program for calling variants in high-throughput sequencing data, 
	particularly paired tumor-normal samples, in 5 distinct steps.
"""

__author__ = "Oliver"
__version__ = "$Revision: 2.0 $"
__date__ = "$Date: 07-2015 $"

import argparse
import sys 
import re
import subprocess
import os
from distutils import spawn

# -------------------------------------

def main():
	"""Main block"""

	# get arguments
	args = get_arg()

	# check for errors, check dependencies
	check_error(args)

	if "1" in args.steps:
		mystep1 = step1(args, "1")
		mystep1.run()

	if "2" in args.steps:
		mystep2 = step1(args, "2")
		mystep2.run()

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
	parser.add_argument("--region","-r",						help="the genomic region to run SAVI on (default: full range) (example: chr1 or chr1:1-50000000)")
	# parser.add_argument("--index",						help="the integer index of the region you want to run SAVI on. By default, 1 refers to the first chromosome in range.txt, 2 to the second, and so on. If you use the partition flag, 1 refers to the first partition, 2 to the second and so on.  (default: 0, which corresponds to the full range)")
	# parser.add_argument("--partition",						help="Number of bases into which to partition the genome.  Use this flag, only if you want to break your genome into regions smaller than the chr length.  If you use this flag, you must also specify an index. For example, \"--partition 50000000 --index 1\" would refer to chr1:1-50000000 for hg19 (default: not used)")
	parser.add_argument("--names",							help="sample names in a comma-delimited list, in the corresponding order of your bam files (default: names are numerical indicies)")
	parser.add_argument("--compsamp","-c",						help="comma-delimited list of colon-delimited indices of samples to compare with savi (default: everything compared to sample 1) (example: 2:1 would compare the second bam file to the first) (example: 2:1,3:1,3:2 would compare the second to the first, the third to the first, and the third to the second)")
	parser.add_argument("--steps",				default="1245",		help="steps to run (default: 1,2,4,5 (i.e., all except prior generation))")
	parser.add_argument("--ann-genome",			default="hg19",		help="name of the SnpEff genome with which to annotate (default: hg19)")
	parser.add_argument("--memory",		type=int,	default=6,		help="the memory for the (SnpEff) Java virtual machine in gigabytes (default: 6)")
	parser.add_argument("--scripts",			default=software,	help="location of scripts dir (directory where this script resides - use this option only if qsub-ing with the Oracle Grid Engine)")
	parser.add_argument("--mindepth",	type=int,	default=6,		help="the min tot read depth required in at least one sample - positions without this wont appear in pileup file (default: 10). Where the filtering occurs: samtools mpileup post-processing")
	parser.add_argument("--minad",		type=int,	default=2,		help="the min alt depth (AD) in at least one sample to output variant (default: 2). Where the filtering occurs: samtools mpileup post-processing")
	parser.add_argument("--mapqual",	type=int,	default=10,		help="skip alignments with mapQ less than this (default: 10). Where the filtering occurs: samtools mpileup")
	parser.add_argument("--s1adpp",		type=int,	default=3,		help="for filtered report, require the sample 1 (normal) alt depth per position to be less than this (default: 3) (note: this is NOT sample1 alt depth of the given alt but, rather, at the given position). Where the filtering occurs: generating report.coding.somatic")
	parser.add_argument("--minallelefreq",		type=int,	default=4,	help="Sgt1MAXFREQ (the allele frequency of any sample not including the first one, assumed to be normal) is greater than this (default: 4) Where the filtering occurs: generating the PD.report file.")
	parser.add_argument("--ann-vcf",						help="comma-delimited list of vcfs with which to provide additional annotation (default: none). Where it's used: SnpSift")
	parser.add_argument("--buildprior",			default=software + "/bin/prior_unif01",		help="starting input prior when building the prior if step 3 is invoked (default: bin/prior_unif01)")
	parser.add_argument("--prior",				default=software + "/bin/prior_diploid01",	help="prior to use if step 3 is not run (default: bin/prior_diploid01)")
	parser.add_argument("--presence",			default="1e-6",		help="the SAVI presence posterior (default: 1e-6). Where it's used: step 4")
	parser.add_argument("--conf",				default="1e-5",		help="the SAVI conf (default: 1e-5). Where it's used: step 4")
	parser.add_argument("--precision",	type=int,	default=0,		help="the SAVI precision (default: 0). Where it's used: step 4")
	parser.add_argument("--rnabams",						help="comma-delimited list of rna bam files (.bai indices should be present)")
	# parser.add_argument("--nofilter",	action="store_true",			help="do not use SAVI comparison filter (default: off) (you should use this option if NOT doing comparisons - e.g., if you only have one sample)")
	parser.add_argument("--nopv4",		action="store_true",			help="do not run bcftools to compute PV4 (default: off). Where it's used: step 5")
	parser.add_argument("--noindeldepth",	action="store_true",			help="do not include include indel reads in total depth count (SDP) (default: off) (note: asteriks in mpileup are included irrespective of this flag). Where it's used: step 2")
	parser.add_argument("--rdplusad",	action="store_true",			help="use reference-agreeing reads plus alternate-calling reads (RD+AD) rather than total depth (SDP) as input to savi (default: off). Where it's used: step 2")
	# parser.add_argument("--hybrid",	action="store_true",			help="as input to savi, use reference-agreeing reads plus alternate-calling reads (RD+AD) for first sample (normal) and SDP for other samples (default: off) (note: this flag changes read depths on positions where there are multiallelic variants)")
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

	# add key-value pairs to the args dict
	vars(args)['numsamp'] = numsamp 
	vars(args)['cwd'] = cwdir
	# if user didn't define compsamp, set it
	if not args.compsamp: vars(args)['compsamp'] = generate_compsamp(numsamp)

	# need to fix paths if users passed software path manually

	# print args
	print(args)

	return args

# -------------------------------------

def check_error(args):
	"""Check for errors, check dependencies """

	# need to check versions
	# need to print times

	# required programs
	required_programs = ("java", "python", "samtools", "snpEff.jar", "bgzip", "tabix", "vcffilter", "bcftools")
	required_savi_programs = ("pileup2multiallele_vcf", "add_to_info", "make_qvt", "savi_poster", "savi_conf", "savi_comp", "savi_poster_accum", "savi_poster_merge", "savi_poster_prt", "savi_unif_prior", "savi_txt2prior")

	# check for required programs 
	for j in required_savi_programs:
		if ( not os.path.isfile(args.scripts + "/bin/" + j) ):
			print("Can't find " + j + ". Please make sure it's in " + software + "/bin")
			sys.exit(1)

	for j in required_programs:
		if (not spawn.find_executable(j)):
			print("Can't find " + j + ". Please add it to your PATH")
			sys.exit(1)

# -------------------------------------

def generate_compsamp(numsamp):
	"""Generate sample comparison string"""

	if (numsamp == 1):
		return '1'
	elif (numsamp > 1):
		return ",".join(["%s:1" % k for k in range(2,numsamp+1)]) 
	else:
		return ''

# -------------------------------------

class step():
	"""A parent step class from which step children inherit"""
	def __init__ (self, args, step_index):
		self.args = args 
		self.step_index = step_index 

	def run(self):
		print('[STEP ' + self.step_index + ']')

# -------------------------------------

class step1(step):
	"""A step1 object whose run method converts bams to mpileup"""

# -------------------------------------

class step2(step):
	"""A step2 object whose run method converts mpileup to vcf"""

# -------------------------------------

class step3(step):
	"""A step3 object whose run method computes prior"""

# -------------------------------------

class step4(step):
	"""A step4 object whose run method runs savi proper"""

# -------------------------------------

class step5(step):
	"""A step5 object whose run method annotates with SnpEff and generates report"""

# -------------------------------------

if __name__ == "__main__":

	main()	
