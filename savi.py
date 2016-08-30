#!/usr/bin/env python

"""
    SAVI
    ~~~~~~
    SAVI is a program for calling variants in high-throughput sequencing data, 
    particularly paired tumor-normal samples, in 5 distinct steps.
"""

__author__ = 'Oliver'
__version__ = 'Revision: 3.0'
__date__ = 'Date: 07-2015'

import argparse, sys, re, subprocess, os, time
from distutils import spawn
# import vcf
from scipy import stats

# -------------------------------------

def main():
    """Main block"""

    # get arguments dictionary
    (args, parser) = get_arg()

    # check for errors, check dependencies
    if not args.noerror:
        check_error(args, parser)

    # make a Cleaner object (keeps a list of files to delete)
    cleaner = Cleaner()

    # if user requests, run the appropriate step
    for i in args.steps:
        # get the class from the classes in the current script
        myclass = globals()['Step' + i]
        # instantiate a step obj for the class
        mystep = myclass(args, i)
        # execute the run method
        mystep.run()
        # update set with list of files to be deleted at the very end
        cleaner.add_step_junk(mystep)

    # clean up tmp files
    if not args.noclean:
        cleaner.cleanup(args.verbose)

# -------------------------------------

def get_arg():
    """Get Arguments"""

    prog_description = '''SAVI is a program for calling variants in high-throughput sequencing data, particularly paired tumor-normal samples. In bioinformatics, \"calling variants\" can mean two different things: (1) simply enumerating all differences between some sequence data and a reference; and (2) determining which of those differences are significant and not likely to be error. SAVI does the later, while a program like Samtools mpileup will do the former. In practice, SAVI is a way of sifting through large amounts of data to pull out the significant mutations using a Bayesian probability model. A common use case is identifying deleterious mutations in cancer, given normal and tumor sequence data---an often-encountered problem in bioinformatics. The output of SAVI is a list of candidate genomic alterations each annotated with a probability of how likely it is to be real. SAVI works with standard bioinformatic file formats. For a complete description of how to use this software, including dependencies and usage examples, see https://github.com/RabadanLab/SAVI'''

    # directory where this script resides
    software = os.path.dirname(os.path.realpath(__file__))
    # cwd
    cwdir = os.getcwd()

    # http://docs.python.org/2/howto/argparse.html
    parser = argparse.ArgumentParser(description=prog_description)

    parser.add_argument('--bams','-b', help='comma-delimited list of bam files; list the normal sample first, as in: normal.bam,tumor.bam (.bai indices should be present)')
    parser.add_argument('--ref', help='reference fasta file with faidx index in the same directory')
    parser.add_argument('--outputdir','-o', default = cwdir, help='the output directory (default: cwd)')
    parser.add_argument('--region','-r', default = '', help='the genomic region to run SAVI on (default: full range) (example: chr1 or chr1:1-50000000)')
    parser.add_argument('--names', help='sample names in a comma-delimited list, in the corresponding order of your bam files (default: names are numerical indicies)')
    parser.add_argument('--compsamp','-c', help='comma-delimited list of colon-delimited indices of samples to compare with savi (default: everything compared to sample 1) (example: 2:1 would compare the second bam file to the first) (example: 2:1,3:1,3:2 would compare the second to the first, the third to the first, and the third to the second)')
    parser.add_argument('--steps', default='1245', help='steps to run (default: 1,2,4,5 (i.e., all except prior generation))')
    parser.add_argument('--ann', default='hg19', help='name of the SnpEff genome with which to annotate (default: hg19)')
    parser.add_argument('--memory', default='4', help='the memory for the (SnpEff) Java virtual machine in gigabytes (default: 6)')
    parser.add_argument('--scripts', default=software, help='location of scripts dir (directory where this script resides - use this option only if qsub-ing with the Oracle Grid Engine)')
    parser.add_argument('--mindepth', type=int, default=10, help='the min tot read depth required in at least one sample - positions without this wont appear in pileup file (default: 10). Where the filtering occurs: samtools mpileup post-processing')
    parser.add_argument('--minad', type=int, default=2, help='the min alt depth (AD) in at least one sample to output variant (default: 2). Where the filtering occurs: samtools mpileup post-processing')
    parser.add_argument('--mapqual', type=int, default=10, help='skip alignments with mapQ less than this (default: 10). Where the filtering occurs: samtools mpileup')
    parser.add_argument('--maxdepth', type=int, default=100000, help='max per-BAM depth option for samtools mpileup (default: 100000)')
    parser.add_argument('--s1adpp', default='3', help='for filtered report, require the sample 1 (normal) alt depth per position to be less than this (default: 3) (note: this is NOT sample1 alt depth of the given alt but, rather, at the given position). Where the filtering occurs: generating report.coding.somatic')
    parser.add_argument('--minallelefreq', default='4', help='Sgt1MAXFREQ (the allele frequency of any sample not including the first one, assumed to be normal) is greater than this (default: 4) Where the filtering occurs: generating the PD.report file.')
    parser.add_argument('--annvcf', help='comma-delimited list of vcfs with which to provide additional annotation (default: none). Where it\'s used: SnpSift')
    parser.add_argument('--buildprior', default=software + '/bin/prior_unif01', help='starting input prior when building the prior if step 3 is invoked (default: bin/prior_unif01)')
    parser.add_argument('--prior', default=software + '/bin/prior_diploid01', help='prior to use if step 3 is not run (default: bin/prior_diploid01)')
    parser.add_argument('--prioriterations', default='10', help='the number of iterations for the prior build, if step 3 is run (default: 10)')
    parser.add_argument('--presence', default='1e-6', help='the SAVI presence posterior (default: 1e-6). Where it\'s used: step 4')
    parser.add_argument('--conf', default='1e-5', help='the SAVI conf (default: 1e-5). Where it\'s used: step 4')
    parser.add_argument('--precision', type=int, default=0, help='the SAVI precision (default: 0). Where it\'s used: step 4')
    parser.add_argument('--noindeldepth', action='store_true', help='do not include include indel reads in total depth count (SDP) (default: off) (note: asteriks in mpileup are included irrespective of this flag). Where it\'s used: step 2')
    parser.add_argument('--rdplusad', action='store_true', help='use reference-agreeing reads plus alternate-calling reads (RD+AD) rather than total depth (SDP) as input to savi (default: off). Where it\'s used: step 2')
    parser.add_argument('--index', default='0', help='an index used in naming of output files (default: 0)')
    parser.add_argument('--rnabams', help='comma-delimited list of rna bam files against which to cross-check your variants (.bai indices should be present)')
    parser.add_argument('--vcf', help='(for running step 5 as a stand-alone only) vcf to use as input')
    parser.add_argument('--noncoding', action='store_true', help='use snpEff to find all transcripts, not just only protein transcripts (default: off). Where it\'s used: step 5')
    parser.add_argument('--noclean', action='store_true', help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--noerror', action='store_true', help='do not check for errors (default: off)')
    parser.add_argument('--verbose','-v', action='store_true', help='echo commands (default: off)')
    parser.add_argument('--superverbose', action='store_true', help='print output of the programs called in each step (default: off)')

    args = parser.parse_args()

    # automatically modify the argument dict:

    # get the number of samples
    numsamp = 0
    if args.bams:
        numsamp = len(args.bams.split(','))

    # if user passes scripts dir, change paths of priors (default and build)
    if args.scripts:
        # change default if user didn't supply his own input
        if args.buildprior == software + '/bin/prior_unif01':
            args.buildprior = args.scripts + '/bin/prior_unif01'
        if args.prior == software + '/bin/prior_diploid01':
            args.prior = args.scripts + '/bin/prior_diploid01'

    # add key-value pairs to the args dict
    vars(args)['numsamp'] = numsamp 
    vars(args)['cwd'] = cwdir
    vars(args)['bin'] = args.scripts + '/bin'
    # if user didn't define compsamp, set it
    if not args.compsamp: vars(args)['compsamp'] = generate_compsamp(numsamp)
    # generate default prior string if user not running step3
    if not '3' in args.steps:
        vars(args)['priorstring'] = generate_priorstr(args.compsamp, args.prior)
    # define a directory for reports
    vars(args)['reportdir'] = args.outputdir + '/report'

    # print args
    print(args)
    print

    return args, parser

# -------------------------------------

def check_error(args, parser):
    """Check for errors, check dependencies"""

    # if no args, print help
    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit(0)

    # check for existence of output dir
    check_dirpath(args.outputdir)

    # unpaired?
    if args.numsamp == 1:
        print('[ALERT] running savi with an unpaired sample - this is not optimal')

    # check step string only goes from 1 to 5 and is all digits 
    if not args.steps.isdigit():
        print('[ERROR] invalid step string')
        sys.exit(1)
    if set(list(args.steps)).intersection(set(list('67890'))):
        print('[ERROR] invalid step string (only steps 1 thro 5 exist)')
        sys.exit(1)

    # define required programs
    required_programs = ('java', 'python', 'samtools', 'snpEff.jar', 'SnpSift.jar', 'bgzip', 'tabix', 'vcffilter')
    required_savi_programs = ('pileup2multiallele_vcf', 'add_to_info', 'make_qvt', 'savi_poster', 'savi_conf', 'savi_comp', 'savi_poster_accum', 'savi_poster_merge', 'savi_poster_prt', 'savi_unif_prior', 'savi_txt2prior', 'filter_pileup')

    # check for required programs 
    for j in required_savi_programs:
        check_path(args.bin + '/' + j)

    for j in required_programs:
        if (not spawn.find_executable(j)):
            print('[ERROR] Can\'t find ' + j + '. Please add it to your PATH')
            sys.exit(1)

    # check for correct versions 
    effpath = spawn.find_executable('snpEff.jar')
    for myprogram, version, mycmd in zip(['samtools', 'SnpEff'], ['1.2', '4.1'], ['samtools --version', 'java -Xmx1G -jar ' + effpath + ' -version']):
        myout = run_cmd(mycmd, 0, 1)
        # print(myout)
        if 'insufficient memory' in myout:
            print('[ERROR] insufficient memory for java')
            sys.exit(1)
        elif not version in myout.split('\n')[0]:
            print('[ERROR] Must use ' + myprogram + ' verion: ' + version)
            sys.exit(1)

    # check for existence of bam files
    if args.bams:
        for j in args.bams.split(','):
            check_path(j)
            check_path(j + '.bai')

    # check for existence of reference
    if args.ref:
        check_path(args.ref)
        check_path(args.ref + '.fai')

    # check for existence of annotating vcf
    if args.annvcf:
        for i in args.annvcf.split(','):
            check_path(i)

    # if Step1 called, but no bam files or ref supplied, throw error
    if '1' in args.steps and ((not args.bams) or (not args.ref)):
        print('[ERROR] Need to supply .bam files and reference file for Step 1')
        sys.exit(1)

    # check that region is in ref
    if args.region:
        with open(args.ref + '.fai', 'r') as f:
            contents = f.read()
            # the region supplied by the user should be in the first column of the ref.fa.fai file
            # region can look like this, e.g., chr1:120000-130000 (just check the part before the colon)
            if not (args.region.split(':')[0] in [x.split()[0] for x in contents.split('\n') if x]):
                print('[ERROR] Region ' + args.region.split(':')[0] + ' not in ' + args.ref + '.fai')
                sys.exit(1)
            
# -------------------------------------

def check_path(myfile):
    """Check for the existence of a file"""

    # expanduser handle tilde
    if (not os.path.isfile(os.path.expanduser(myfile))):
        print('[ERROR] Can\'t find the file ' + myfile)
        sys.exit(1)

# -------------------------------------

def check_dirpath(myfile):
    """Check for the existence of a directory"""

    # expanduser handle tilde
    if (not os.path.isdir(os.path.expanduser(myfile))):
        print('[ERROR] Can\'t find the folder ' + myfile)
        sys.exit(1)

# -------------------------------------

def generate_compsamp(numsamp):
    """Generate sample comparison string, which determines which samples are compared to which"""

    # by default, everything is compared to sample 1, which is assumed to be 'normal' (i.e., not tumor)
    if numsamp == 1:
        return '1'
    elif numsamp > 1:
        return ','.join(['%s:1' % k for k in range(2,numsamp+1)])
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
    # if compsamp='2:1', it should look like this - 1:${outputdir}/savi/prior1/prior,2:${outputdir}/savi/prior2/prior
    # if not step3, use diploid prior by default
    if not compsamp:
        return '1:' + priorpath
    else:
        return ','.join(['%s:%s' % (k, priorpath) for k in get_uniq_samples(compsamp)])

# -------------------------------------

def run_cmd(cmd, bool_verbose, bool_getstdout):
    """Run a system (i.e., shell) command"""

    # if verbose, print command
    if (bool_verbose):
        print('[command] ' + cmd)
        sys.stdout.flush()

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() 
    # print return code
    # print(proc.returncode) 
    # print stdout stderr tuple
    # proc.communicate()

    (stdout, stderr) =  proc.communicate()

    # if error, print it
    if stderr:
        print('[stderror] ' + stderr),

    # return stdout
    if (bool_getstdout): 
        return stdout.rstrip()
    else:
        # note: this must return a str
        return '0'

# -------------------------------------

def run_cmd_python_pipe(cmd, bool_verbose, filter_function=None, output_file=None):
    """Run a system (i.e., shell) command and filter it with python as if we piped it into python"""

    # if verbose, print command
    if (bool_verbose):
        print('[command] ' + cmd)

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # http://stackoverflow.com/questions/2804543/read-subprocess-stdout-line-by-line

    with open(output_file, 'w') as f:
        while True:
            line = proc.stdout.readline()
            if line:
                myline = filter_function(line.rstrip())
                if myline:
                    f.write(myline + '\n')
            else:
                break

# -------------------------------------

def check_file_exists_and_nonzero(myfile):
    """Check for the existence and nonzero-ness of a file"""

    # loop through comma-delimited list of files
    for i in myfile.split(','):
        if os.path.isfile(os.path.expanduser(i)):
            if os.path.getsize(os.path.expanduser(i)) == 0:
                print(i + ' is empty. Exiting')
                sys.exit(1)
        else:
            print('Can\'t find ' + i + '. Exiting.')
            sys.exit(1)

# -------------------------------------

def mytimer(myfunc):
    """Decorator for timing a Step's run function"""
    # http://stackoverflow.com/questions/5478351/python-time-measure-function

    def mynewfunc(*args, **kwargs):
        startTime = time.time()
        myfunc(*args, **kwargs)
        print('[step delta t] {} sec'.format(int(time.time() - startTime)))

    return mynewfunc

# -------------------------------------

# not used
def filter_pileup(args, myline):
    """Filter line of pileup"""
    # reduce pileup by only printing variants and only printing lines with enough depth

    # number of samples
    num = int(args.numsamp)
    # tot read cutoff
    cutoff = int(args.mindepth)
    # alt read depth
    minad = int(args.minad)

    # if sample 1 (normal) doesn't meet cutoff, break
    if int(myline.split()[3]) < cutoff:
        return None

    # if only one sample, print
    if num == 1:
        # print(myline)
        return myline
    # if multiple samples
    else:
        # loop thro samples, discounting the normal sample (sample 1)
        for i in range(1,num):
            # only consider if sample depth greater than or eq to cutoff
            if int(myline.split()[3+i*3]) >= cutoff:
                # if make prior step, generate all positions in the mpileup, variants and non-variants alike
                if '3' in args.steps:
                    return myline
                    break
                else:
                    # get read string
                    readstr = myline.split()[4+i*3]
                    # count ACTGN mismatches
                    # if sufficient number of mismatches, flag line to print
                    # but eliminate read starts with qual scores ACTG (e.g., stuff like ^A ^C etc)
                    if len(re.findall(r'[AaGgCcTtNn]', re.sub(r'\^[ACTGNactgn\+\-]', '', readstr))) >= minad:
                        return myline
                        break

    return None

# -------------------------------------

class Cleaner():
    """A class whose object keeps a set of files to delete from the Steps"""

    def __init__ (self):
        """Initialize Cleaner object"""

        # set of intermediate files (to be cleaned)
        self.junk=set()

    def add_step_junk(self, myStep):
        """Update self.junk with the list of junk files from the step object"""

        # update set with list
        self.junk.update(myStep.intermediates)

    def cleanup(self, bool_verbose):
        """Clean up tmp and intermediate files"""

        if bool_verbose:
            print('[clean-up] Remove files: '),
            print(self.junk)

        # rm -f
        for i in self.junk:
            if os.path.isfile(os.path.expanduser(i)):
                os.remove(i)

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

        # extra carriage return first for all steps other than 1
        if self.args.verbose and self.step_index != '1':
            print

        # print step and description
        print(('[STEP {self.step_index}] {self.description}').format(self=self))

# -------------------------------------

class Step1(Step):
    """A step1 object whose run method converts bams to mpileup"""

    def __init__ (self, args, step_index):
        """Step1 constructor"""

        # parent's constructor
        Step.__init__(self, args, step_index)

        # define description as well as input and output attributes
        self.set_descrip('Convert bams to mpileup and, if not building custom prior, filter for variants only')
        self.set_input(self.args.bams)
        self.set_output(('{self.args.outputdir}/tmp_mpile.{self.args.index}.txt').format(self=self))
        # for rna bams:
        self.output2 = '{self.args.outputdir}/tmp_mpile.{suffix}.txt'.format(self=self, suffix=int(self.args.index) + 1)

    # override parent's run method, decorate with mytimer to benchmark time
    @mytimer
    def run(self):
        """The run method calls shell(system) commands to do step1 - namely, it calls samtools mpileup to convert bam to mpileup"""

        # run parent's run method
        super(Step1, self).run()

        # define pileup file
        # note the -d flag: 'max per-BAM depth to avoid excessive memory usage [250].'
        # The default of 250 is too low, so we jack it up
        pileupflag='-A -B -q {self.args.mapqual} -d {self.args.maxdepth} -L {self.args.maxdepth} -f {self.args.ref}'.format(self=self)

        # define command to run

        # filer pileup such that: normal has depth > cutoff
        # and at least one of the tumor samples has depth > cutoff
        # and we want only variants
        # (this implemented in Python in the function filter_pileup above, but can't figure out how to make it run fast so use small binary)
        filtercmd = '{self.args.bin}/filter_pileup --minaltdepth {self.args.minad} --mindepth {self.args.mindepth} --normaldepth'.format(self=self)

        # if make prior step, filter to generate all positions in the mpileup, variants and non-variants alike
        if '3' in self.args.steps:
            filtercmd = '{self.args.bin}/filter_pileup --mindepth {self.args.mindepth} --buildprior'.format(self=self)

        # define region flag
        regionflag = ''
        if self.args.region:
            regionflag='-r ' + self.args.region

        # command: samtools plus filter pileup by only printing variants and only printing lines with enough depth
        mycmd = 'samtools mpileup {} {} {} | {} > {}'.format(pileupflag, regionflag, self.args.bams.replace(',', ' '), filtercmd, self.output)
        # run it
        run_cmd(mycmd, self.args.verbose, 1)

        # check if output file nonzero size
        check_file_exists_and_nonzero(self.output)

        # the output of this step is needed for the next step, but ultimately we want to delete it
        self.intermediates.append(self.output)

        # repeat with slightly different filtercmd for RNA
        if self.args.rnabams:
            # define filtercmd
            filtercmd = '{self.args.bin}/filter_pileup --minaltdepth {self.args.minad} --mindepth {self.args.mindepth}'.format(self=self)
            mycmd = 'samtools mpileup {} {} {} | {} > {}'.format(pileupflag, regionflag, self.args.rnabams.replace(',', ' '), filtercmd, self.output2)
            run_cmd(mycmd, self.args.verbose, 1)
            self.intermediates.append(self.output2)

# -------------------------------------

class Step2(Step):
    """A step2 object whose run method converts mpileup to vcf"""

    def __init__ (self, args, step_index):
        """Step2 constructor"""

        # parent's constructor
        Step.__init__(self, args, step_index)

        # define description as well as input and output attributes
        self.set_descrip('Convert mpileup to vcf')
        self.set_input(self.args.outputdir + '/' + 'tmp_mpile.' + self.args.index + '.txt')
        self.set_output(self.args.outputdir + '/' + self.args.index + '.vcf.bgz')
        # for rna bams:
        self.input2 = self.args.outputdir + '/' + 'tmp_mpile.' + str(int(self.args.index) + 1) + '.txt'
        self.output2 = self.args.outputdir + '/' + str(int(self.args.index) + 1) + '.vcf'

    # override parent's run method
    @mytimer
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

        # RNA
        if self.args.rnabams:
            mycmd = "cat {} | {}/pileup2multiallele_vcf {} > {}".format(self.input2, self.args.bin, oscanflag, self.output2)
            run_cmd(mycmd, self.args.verbose, 1)

        # if empty
        # if [ $( zcat ${outputdir}/${SGE_TASK_ID}.vcf.bgz | sed '/^#/d' | head | wc -l ) == 0 ]; then
        #     echo "[HALT SCRIPT] vcf file is empty"
        #     exit 0;
        # fi

        # the output of this step is needed for the next step, but ultimately we want to delete it
        for i in [self.output, self.output + ".tbi", self.output2]:
            self.intermediates.append(i)

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

    @mytimer
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

        # vcfs
        self.vcf_freq = self.args.reportdir + "/freqsavi.vcf.bgz"
        self.vcf_final = self.args.reportdir + "/finalsavi.vcf.bgz"

    # override parent's run method
    @mytimer
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
            os.rename(self.vcf_freq, self.vcf_final)
            os.rename(self.vcf_freq + '.tbi', self.vcf_final + '.tbi')

        mycmd = "gunzip -S .bgz " + self.vcf_final
        myout = run_cmd(mycmd, self.args.verbose, 1)

        # the output of this step is needed for the next step, but ultimately we want to delete it
        # fix this hack later
        for i in [self.vcf_freq, self.vcf_freq.replace(".bgz",""), self.vcf_freq + ".tbi", self.vcf_final, self.vcf_final.replace(".bgz",""), self.vcf_final + ".tbi"]:
            self.intermediates.append(i)

# -------------------------------------

class Step5(Step):
    """A step5 object whose run method annotates with SnpEff and generates report"""

    def __init__ (self, args, step_index):
        """Step5 constructor"""

        # parent's constructor
        Step.__init__(self, args, step_index)

        # define description as well as input and output attributes
        self.set_descrip("Add annotation to the vcf")
        if self.args.vcf:
            self.set_input(self.args.vcf)
        else:
            self.set_input(self.args.reportdir + "/finalsavi.vcf")

        # set output prefix, output reports
        self.outprefix = self.args.reportdir + "/tmp_" + self.args.index

        # full vcf report 
        self.vcf_all = self.args.reportdir + "/report.all.vcf"
        # vcf report, filtered for coding regions 
        self.vcf_coding = self.args.reportdir + "/report.coding.vcf"
        # vcf report, filtered for coding regions + somatic variants
        self.vcf_somatic = self.args.reportdir + "/report.coding.somatic.vcf"
        # vcf report, add PV4
        self.vcf_somatic_pv4 = self.args.reportdir + "/report.coding.somatic.pv4.vcf"
        # vcf report, filtered for PD_L > 0
        self.vcf_PD = self.args.reportdir + "/report.coding.PDfilter.vcf"
        # vcf report, filtered for PD_U < 0
        self.vcf_PD_rev = self.args.reportdir + "/report.coding.rev-PDfilter.vcf"

    def runEff(self):
        """Run SnpEff and SnpSift"""

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
        
        # add file to list of intermediates
        self.intermediates.append(self.output)

        # define SnpSift commands to run

        # loop through all annotating vcfs
        if self.args.annvcf:
            for j,i in enumerate(self.args.annvcf.split(",")):
                # set output file
                self.output = self.outprefix + ".eff_" + str(j+1) + ".all.vcf"
                # add file to list of intermediates
                self.intermediates.append(self.output)

                mycmd = "java -Xmx" + self.args.memory + "G" + \
                    " -jar " + siftpath + \
                    " annotate " + i + \
                    " -noLog " + \
                    self.outprefix + ".eff_" + str(j) + ".all.vcf > " + self.output

                myout = run_cmd(mycmd, self.args.verbose, 1)
                if (self.args.superverbose): print(myout)

        # rename final iteration of vcf to report
        if (self.args.verbose): print("[comment] rename " + self.output + " to " + self.vcf_all)
        os.rename(self.output, self.vcf_all)

    def generate_PD_str(self, mycompsamp):
        """generate PD filter string"""

        # initialize output strings
        outstr=""
        outstr_rev=""

        # get comparisons only
        compsamp_list = [i.replace(":", "") for i in mycompsamp.split(',') if ':' in i]

        # generate filter string for vcffilter
        if compsamp_list:
            for j,k in enumerate(compsamp_list):
                if (j < len(compsamp_list) - 1):
                    outstr += "PD" + k + "_L > 0 | "
                    outstr_rev += "PD" + k + "_U < 0 | "
                else:
                    outstr += "PD" + k + "_L > 0"
                    outstr_rev += "PD" + k + "_U < 0"

        return(outstr, outstr_rev)
        
    def filterNs(self, report_in, report_out):
        """Filter out N chars in vcf ref and alt"""

        if (self.args.verbose): print("[comment] filter Ns from " + report_in + ", produce: " + report_out)

        # make sure neither ref nor alt column has an 'N'
        with open(report_out, 'w') as g:
            with open(report_in, 'r') as f:
                for line in f:
                    if (line.startswith("#") or ((not 'N' in line.split()[3]) and (not 'N' in line.split()[4]))):
                        g.write(line)

    def filterCoding(self, report_in, report_out):
        """Filter for the coding region - 
        get useful variants and variants in the coding region. The idea is to capture these SnpEff features:
        inframe_deletion *inframe_insertion frameshift_variant* initiator_codon_variant missense_variant* synonymous_variant*
        splice_acceptor* splice_donor* splice_region* start_lost* start_gained* stop_gained* stop_lost* stop_retained_variant*;
        Also, add strand bias p-value"""

        if (self.args.verbose): print("[comment] filter " + report_in + " for coding region, produce: " + report_out)

        # filter for coding region
        # desired features
        featurelist = ["inframe", "frameshift", "synonymous_variant", "missense_variant", "splice_acceptor", "splice_donor", "splice_region", "start", "stop"]
        with open(report_out, "w") as g:
            with open(report_in, 'r') as f:
                for line in f:
                    # header line
                    if line.startswith("#"): 
                        # write to file
                        g.write(line)
                        # this is a hack to add the strand bias into the header
                        if line.startswith('##FORMAT=<ID=ADR,'):
                            g.write('##FORMAT=<ID=SB,Number=1,Type=Float,Description="P-value for strand bias (exact test)">\n')
                    # if non-header line with desired features 
                    elif any(k in line.split()[7] for k in featurelist):

                        try:
                            # print line, but add strand bias p-value

                            # get indices of forward and reverse read counts
                            irdf = line.split()[8].split(":").index('RDF')
                            irdr = line.split()[8].split(":").index('RDR')
                            iadf = line.split()[8].split(":").index('ADF')
                            iadr = line.split()[8].split(":").index('ADR')
                        
                            # variable for format fields, so can add strand bias
                            format_fields = ""

                            # calculate strand bias
                            for i in range(9,len(line.split())):
                                oddsratio, pvalue = stats.fisher_exact([[line.split()[i].split(':')[irdf], line.split()[i].split(':')[irdr]], [line.split()[i].split(':')[iadf], line.split()[i].split(':')[iadr]]])
                                format_fields += line.split()[i] + ":" + str(round(pvalue,6)) + "\t"

                            # write to file
                            g.write("\t".join(line.split()[0:9]) + ":SB\t" + format_fields.strip() + "\n")
                        except:
                            g.write(line)

    def filterCodingSomatic(self, report_in, report_out):
        """Filter for the coding + somatic variants"""

        if (self.args.verbose): print("[comment] filter " + report_in + " for somatic variants and discard meganormals and common dbSNps, produce: " + report_out)

        # define command
        # filtering: S1 AD PP < cutoff so discard variants with normal alt depth above a certain threshhold
        # discard variants in meganormals
        # discard variants in COMMON dbSNP
        # This is awk-ward (pun intended!) but use for now
        mycmd = 'vcffilter -f "S1ADPP < ' + self.args.s1adpp + '" ' + report_in + " | awk '$0 ~ /^#/ || ( $0 !~ /MEGANORMAL_ID/ && $0 !~ /Meganormal186GBM/ && $0 !~ /COMMON\=1/)' > " + report_out
        myout = run_cmd(mycmd, self.args.verbose, 1)

    def filterCodingPD(self, report_in, report_out, report_out_rev):
        """Filter for the PD lower bound greater than zero, upper bound less than zero"""

        if (self.args.verbose): print("[comment] filter " + report_in + " for PD lower bound > 0 and PD upper bound < 0, produce: " + report_out + ", " + report_out_rev)

        # get filter strings
        filterstring, filterstring_reverse = self.generate_PD_str(self.args.compsamp)

        # if filterstring defined
        if filterstring:
            # add filter for PD00 if at least 3 samples
            if self.args.numsamp > 2: 
                filterstring = "( Sgt1MAXFREQ > " + self.args.minallelefreq + " ) & ( ( P1_L > 0 & P1_L < 51 & P1_U > 49 ) | ( P1_L = 0 ) ) & ( PD00_U < 0 | " + filterstring + " )"
                filterstring_reverse = "PD00_L > 0 | " + filterstring_reverse
            else:
                filterstring = "( P2_F > " + self.args.minallelefreq + " ) & ( ( P1_L > 0 & P1_L < 51 & P1_U > 49 ) | ( P1_L = 0 ) ) & ( " + filterstring + " )"
                filterstring_reverse = "( P1_L > 0 & P1_L < 51 & P1_U > 49 ) & ( " + filterstring_reverse + " )"

            # TO DO: explain the logic of these filters

            # define command
            mycmd = 'vcffilter -f "' + filterstring + '" ' + report_in + " > " + report_out
            myout = run_cmd(mycmd, self.args.verbose, 1)
            mycmd = 'vcffilter -f "' + filterstring_reverse + '" ' + report_in + " > " + report_out_rev
            myout = run_cmd(mycmd, self.args.verbose, 1)

    def vcf2tsv(self, report_list_in):
        """Turn vcf into human readable, tab delimited report"""

        # if sample names given, use them
        reportoptions = ""

        if self.args.names:
            reportoptions = "--samples " + self.args.names

        for i in report_list_in:
            if os.path.isfile(os.path.expanduser(i)):
                # the sed nonsense is so Excel renders it properly
                mycmd = "cat " + i + " | " + self.args.bin + "/vcf2newreport_for_v4.1.C.py " + reportoptions + " | sed 's|,|, |g' > " + i.replace('.vcf', '.txt')
                myout = run_cmd(mycmd, self.args.verbose, 1)

                # if rna bams present, add depths of rna variants
                if self.args.rnabams:
                    mycmd = self.args.bin + "/join_tsv2vcf.py " + i.replace('.vcf', '.txt') + " " + self.args.outputdir + "/" + str(int(self.args.index) + 1) + ".vcf > " + i.replace('.vcf', '.rna.txt')
                    myout = run_cmd(mycmd, self.args.verbose, 1)

    def addPV4(self, report_in, report_out):
        """Add PV4, defined by samtools as "P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t); 4) tail distance bias (t)" """
        # see http://samtools.sourceforge.net/mpileup.shtml

        # http://pyvcf.readthedocs.org/en/latest/INTRO.html
        # Annoyingly, "in general modifying and writing VCF files is not covered well (yet) in PyVCF" 
        # https://github.com/jamescasbon/PyVCF/issues/82
        # so stop using it and do it by hand ! 

        #vcf_reader = vcf.Reader(open(myvcf, 'r'))
        #for record in vcf_reader:
        #    print record
        #    for j in record.samples:
        #        oddsratio, pvalue = stats.fisher_exact([[j['RDF'], j['RDR']], [j['ADF'], j['ADR']]])
        #        print(pvalue)

        pass

    def runPostEffProcessing(self):
        """Do the Post-SnpEff Processing - Filter vcf and convert to human readable report"""

        # convert vcf to tsv (tab-delimited report)
        mycmd = "cat " + self.vcf_all + " | " + self.args.bin + "/vcf2fullreport.py > " + self.vcf_all.replace('.vcf', '.txt')
        run_cmd(mycmd, self.args.verbose, 1)

        # filter for coding region
        self.filterCoding(self.vcf_all, self.vcf_coding)
        # filter for somatic region
        self.filterCodingSomatic(self.vcf_coding, self.vcf_somatic)
        # TESTING
        self.addPV4(self.vcf_somatic, self.vcf_somatic_pv4)
        # PD filter
        self.filterCodingPD(self.vcf_coding, self.vcf_PD, self.vcf_PD_rev)
        # convert vcfs to human readable tsv reports
        self.vcf2tsv([self.vcf_coding, self.vcf_somatic, self.vcf_PD, self.vcf_PD_rev])

    @mytimer
    def run(self):
        """The run method calls shell(system) commands to do step5 - namely, annotate the vcf"""

        # run parent's run method
        super(Step5, self).run()

        # get rid of Ns before annotation because SnpEff messes this up (see bug report: https://github.com/pcingola/SnpEff/issues/54)
        report_no_N = self.args.reportdir + "/finalsavi.no_N.vcf" 
        self.filterNs(self.input, report_no_N)
        # update input attribute
        self.input = report_no_N
        # add no_N report to list of intermediates
        self.intermediates.append(report_no_N)

        # run SnpEff
        self.runEff()
        # run post-SnpEff filtering, processing
        self.runPostEffProcessing()

# -------------------------------------

if __name__ == "__main__":

    main()
