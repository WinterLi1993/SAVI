#!/bin/bash

### default variables
outputdir=.					# output directory
java_memory=6					# memory of the Java virtual mach in Gigabytes (default: 6)
SGE_TASK_ID=0					# index of region (default: all (0))
myregion=""					# user specified region (default NULL)
mypartition=""					# user specified partition (default NULL)
debug=0						# bool to turn on debugging (default: no)
filter=1					# bool to do SAVI filter (default: yes)
stepstr=12345					# what steps to run (default: 1 thro 5)
software=$( dirname $( readlink -m $0 ) )	# directory where to look for scripts
input_bams=""					# input bam files
num_bams=2					# number of input bams (default: 2)
compsamp="2:1"					# the indices of sample to compare with savi (default: 2:1 - i.e., tumor vs norm)
cutoff=5 					# min read depth cutoff
anngenome="GRCh37.71"				# annotating reference for SnpEff (default: hg19)
annvcf=""					# annotating vcfs in comma-delimited list (default: none)
						# common hg19 annotating vcfs in comma-delimited list
savi_precision=0				# the SAVI precision, from the Hossein fix
helpmessage=$( cat <<EOF
Usage:

$0 --bam [list_of_bam_files] --ref [ref]

Required Arguments:

  --bam		comma-delimited list of bam files (order should be: normal.bam, tumor.bam) (.bai indices should be present)
  --ref		faidx-indexed ref

Options:

  --outputdir	the output directory (default: cwd)
  --region	the genomic region to run SAVI on (default: full range) (example: chr1 or chr1:1-50000000)
  		If you use this flag, the index you supply becomes the file suffix.
  --index	the integer index of the region you want to run SAVI on.
  		By default, 1 refers to the first chromosome in range.txt, 2 to the second, and so on.
		If you use the partition flag, 1 refers to the first partition, 2 to the second and so on.
 		(default: 0, which corresponds to the full range)
  --partition	Number of bases into which to partition the genome.
  		Use this flag, only if you want to break your genome into regions smaller than the chr length.
		If you use this flag, you must also specify an index. For example, 
	        "--partition 50000000 --index 1" would refer to chr1:1-50000000 for hg19
 		(default: not used)
  --compsamp	comma- colon- delimited indices of samples to compare with savi (default: 2:1) (example: 2:1,3:1,3:2)
  --steps	steps to run (default: 1,2,3,4,5,6 (i.e., all except cnv))
  --memory	the memory for the Java virtual machine in gigabytes (default: 6)
  --mindepth	the minimum read depth required in at least one sample - positions without this wont appear in pileup file (default: 5)
  --debug	do not delete tmp intermediate files  (default: off)
  --nofilter	do not use SAVI comparison filter (default: off) (you should use this option if NOT doing comparisons)
  --scripts	location of scripts dir  (directory where this script resides - use this option only if qsub-ing)
  --ann-genome  name of the SnpEff genome with which to annotate (default: GRCh37.71)
  --ann-vcf     comma-delimited list of vcfs with which to provide additional annotation (default: none)
  --prior	input prior (default: savi/unif_prior01)
  --precision	the SAVI precision
  --help	print this message and exit

Example:

$0 --bam normal.bam,tumor.bam --outputdir . --index 1 --ref /my/ref.fasta

Notes:

This scripts assumes that you have a paired normal-tumor dataset which is mapped to human hg19.
If this is not so, you must change the genome reference as well as the genome partition.
This scripts also assumes that Samtools, bgzip, tabix are in your PATH.
Also, this script uses hard-wired paths to VarScan and SnpEff.
You should change these at the top of this script to reflect their paths on your system.

Your reference should be indexed so a .fai file resides in its directory.
Your bam files should be indexed so .bai files reside in their directories.

Steps:

  1 mpileup
  2 mpileup to multiallelic vcf
  3 make prior
  4 run savi
  5 snpeff
  6 cnv and loh (ONLY WORKS FOR PAIRED NORMAL-TUMOR SAMPLES)

EOF
)

# If no arguments, echo help message and quit
if [ $# == 0 ]; then
	echo "$helpmessage"
	exit;
fi

### getopts 
while [ $# -gt 0 ]; do
	if [  "$1" == "-h" -o "$1" == "-help" -o "$1" == "--help" ]; then
		shift; 
		echo "$helpmessage"
		exit;
	elif [  "$1" == "-outputdir" -o "$1" == "--outputdir" ]; then
		shift; 
		outputdir=$1; 
		shift
	elif [  "$1" == "-scripts" -o "$1" == "--scripts" ]; then
		shift; 
		software=$1; 
		shift
	elif [  "$1" == "-ref" -o "$1" == "--ref" ]; then
		shift; 
		ref=$( readlink -m $1 ); 
		shift
		if [ ! -e $ref ]; then
			echo "Cannot find the reference file:" 
			echo ${ref}
			exit 1
		fi
	elif [  "$1" == "-vcf" -o "$1" == "--vcf" ]; then
		shift; 
		vcf=$1; 
		shift
	elif [  "$1" == "-compsamp" -o "$1" == "--compsamp" ]; then
		shift; 
		compsamp=$1; 
		shift
	elif [  "$1" == "-steps" -o "$1" == "--steps" ]; then
		shift; 
		stepstr=$1; 
		shift
	elif [  "$1" == "-memory" -o "$1" == "--memory" ]; then
		shift; 
		java_memory=$1;
		shift
	elif [  "$1" == "-mindepth" -o "$1" == "--mindepth" ]; then
		shift; 
		cutoff=$1;
		shift
	elif [  "$1" == "-ag" -o "$1" == "--ann-genome" ]; then
		shift; 
		anngenome=$1;
		shift
	elif [  "$1" == "-av" -o "$1" == "--ann-vcf" ]; then
		shift; 
		annvcf=$1;
		shift
	elif [  "$1" == "-debug" -o "$1" == "--debug" ]; then
		shift; 
		debug=1; 
	elif [  "$1" == "-nofilter" -o "$1" == "--nofilter" ]; then
		shift; 
		filter=0; 
	elif [  "$1" == "-region" -o "$1" == "--region" ]; then
		shift; 
		myregion=$1; 
		shift
	elif [  "$1" == "-index" -o "$1" == "--index" ]; then
		shift; 
		SGE_TASK_ID=$1; 
		shift
	elif [  "$1" == "-partition" -o "$1" == "--partition" ]; then
		shift; 
		mypartition=$1; 
		shift
	elif [  "$1" == "-prior" -o "$1" == "--prior" ]; then
		shift; 
		inputprior=$1; 
		shift
	elif [  "$1" == "-precision" -o "$1" == "--precision" ]; then
		shift; 
		savi_precision=$1; 
		shift
	elif [  "$1" == "-bam" -o "$1" == "--bam" ]; then
		shift; 
		input_bams=$1; 
		# replace commas w spaces
		input_bams=$( echo $input_bams | sed 's|,| |g' ) 
		num_bams=$( echo $input_bams | wc -w )
		shift
	else	# if unknown argument, just shift
		shift
	fi
done

time1=$( date "+%s" )
inputprior=${software}/savi/unif_prior01	# input prior (default: savi/unif_prior01)

echo "[start]"
echo "[pwd] "`pwd`
echo "[date] "`date`
echo "[index] "$SGE_TASK_ID
echo "[output_dir] "$outputdir
echo "[steps] "$stepstr
echo "[debug] "$debug
echo "[bam files] "$input_bams
echo "[scripts dir] "$software

# check for bams
for i in $input_bams; do 
	if [ ! -e ${i} ]; then 
		echo "Cannot find the bam file:"
		echo "${i}"
		exit 1
	fi; 
done
# check for ref
if [ ! -e ${ref} ]; then 
	echo "Cannot find the reference file:"
	echo "${ref}"
	exit 1
fi
# check for faidx-indexed ref
if [ -e ${ref}.fai ]; then
	echo "[reference] "$ref
else
	echo "This script requires your reference to be faidx-indexed but it cannot find the file:"
	echo "${ref}.fai"
	echo
	echo "Please run:"
	echo "samtools faidx ${ref}"
	exit 1
fi

mkdir -p ${outputdir}

# turn on debugging
set -eux

# if region set 
if [ ! -z "${myregion}" ]; then
	regionflag="-r ${myregion}";
	echo "[region] "$myregion
# not set
else
	if [ $SGE_TASK_ID == 0 ]; then
		regionflag="";
		echo "[region] full range"
	else
		# if partition set
		if [ ! -z "${mypartition}" ]; then
			# create a file of ranges, broken by partitions (the fancy perl stuff is just to sort the chrs by number for number-string combos)
			cat ${ref}.fai | perl -ne 'BEGIN{@myarr=()}{chomp($_); @a=split(/\t/, $_); push(@myarr, $a[0]."\t".$a[1])}END{my @sorted = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_=~/(\d+)/] } @myarr ; foreach $elt (@sorted) {print $elt."\n"}}' | ${software}/util/partition_genome.pl ${mypartition} > ${outputdir}/range.txt
		else
			# create a file of ranges of chrs (not broken by partitions)
			cat ${ref}.fai | perl -ne 'BEGIN{@myarr=()}{chomp($_); @a=split(/\t/, $_); push(@myarr, $a[0]."\t".$a[1])}END{my @sorted = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_=~/(\d+)/] } @myarr ; foreach $elt (@sorted) {print $elt."\n"}}' | cut -f1 > ${outputdir}/range.txt
		fi

		echo "[region] "$( cat ${outputdir}/range.txt | head -${SGE_TASK_ID} | tail -1 )
		regionflag="-r "$( cat ${outputdir}/range.txt | head -${SGE_TASK_ID} | tail -1 )

	fi
fi

### need to remove duplicate reads w Picard
# or else get false positives

if [[ $stepstr == *1* ]]; then
	echo "[STEP1] pileup"
	### mpileup
	# for mpileup - note the -d flag: "max per-BAM depth to avoid excessive memory usage [250]." The default of 250 is way too low, so we jack it up
	# for awk - we want only lines where normal has depth, and at least one of the tumor samples has depth
	# for awk 2 - we want only variants
	# for perl - varscan has a bug that, if the pileup is 0\t\t , as opposed to 0\t*\t*, VarScan fails to parse it properly
	# cutoff is min read depth cutoff
	#samtools mpileup -d 100000 -L 100000 ${regionflag} -f $ref ${input_bams} \
	#    | awk -v num=${num_bams} -v cutoff=${cutoff} '{true=0; for (i = 1; i <= num-1; i++) {if ($(4+i*3) > cutoff) {true=1}}; if (num==1 && $4>cutoff) {print} else if (num>1 && $4>cutoff && true) {print}}' \
	#    | awk -v num=${num_bams} '{true=0; for (i = 0; i < num; i++) {if (match($(5+i*3),/[ATGC\+\-atgc]/)) {true=1}}; if (true) {print}}' | perl -ne '{ if ($_ =~ s/\t\t/\t\*\t\*/g) {print} else {print} }' \
	#    > ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt
	samtools mpileup -d 100000 -L 100000 ${regionflag} -f $ref ${input_bams} \
	    | awk -v num=${num_bams} -v cutoff=${cutoff} '{true=0; for (i = 1; i <= num-1; i++) {if ($(4+i*3) > cutoff) {true=1}}; if (num==1 && $4>cutoff) {print} else if (num>1 && $4>cutoff && true) {print}}' \
	    > ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt

	echo "[date 1] "`date`;

	# if empty
	if [ ! -s ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt ]; then
		echo "[halt script] pileup file is empty"
		exit 0;
	fi
fi

if [[ $stepstr == *2* ]]; then
	echo "[STEP2] pileup2vcf"
	### Oscan
	# cat ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt | ${software}/savi/pileup2multiallele_vcf --header | bgzip > ${outputdir}/${SGE_TASK_ID}.vcf.bgz
	# call everything, variants and non-variants alike
	cat ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt | ${software}/savi/pileup2multiallele_vcf --all --header \
             | tee ${outputdir}/${SGE_TASK_ID}.all.vcf | awk '{ if ($0 ~ /^#/ || toupper($4) != toupper($5)) {print}}' | bgzip > ${outputdir}/${SGE_TASK_ID}.vcf.bgz

	cat ${outputdir}/${SGE_TASK_ID}.all.vcf | bgzip > ${outputdir}/${SGE_TASK_ID}.all.vcf.bgz
	rm ${outputdir}/${SGE_TASK_ID}.all.vcf

	tabix -p vcf ${outputdir}/${SGE_TASK_ID}.all.vcf.bgz
	tabix -p vcf ${outputdir}/${SGE_TASK_ID}.vcf.bgz
	echo "[date 2] "`date`

	# if empty
	if [ $( zcat ${outputdir}/${SGE_TASK_ID}.vcf.bgz | sed '/^#/d' | head | wc -l ) == 0 ]; then
		echo "[halt script] vcf file is empty"
		exit 0;
	fi

fi

if [[ $stepstr == *3* ]]; then
	echo "[STEP3] make_prior"
	### make_prior.py
	# for samples 1 through n
	#for i in $( seq 1 $num_bams ); do
	# for samples in the compsamp string
	for i in $( echo ${compsamp} | sed 's|,|\n|g; s|:|\n|g' | sort -u ); do
		mkdir -p ${outputdir}/savi_${SGE_TASK_ID}/prior${i}
		${software}/make_prior.py --verbose \
		      --name s${i} \
		      --iteration 10 \
		      --prior ${inputprior} \
		      --sampleindex ${i} \
		      --outputdir ${outputdir}/savi_${SGE_TASK_ID}/prior${i} \
		      --input ${outputdir}/${SGE_TASK_ID}.all.vcf.bgz;
	done;

	echo "[date 3] "`date`;
fi

if [[ $stepstr == *4* ]]; then
	echo "[STEP4] run_savi"
	### run_savi.py
	mkdir -p ${outputdir}/savi_${SGE_TASK_ID}/out

	# prior string for priors specified in compsamp
	# if compsamp="2:1", it should look like this - 1:${outputdir}/savi_${SGE_TASK_ID}/prior1/prior,2:${outputdir}/savi_${SGE_TASK_ID}/prior2/prior
	priorstring=$( for i in $( echo ${compsamp} | sed 's|,|\n|g; s|:|\n|g' | sort -u ); do echo -n "${i}:${outputdir}/savi_${SGE_TASK_ID}/prior${i}/prior,"; done | sed 's|,$||'; echo )

	# invert filter variable
	keepfreqfile=$( echo $filter | awk '{x=$1; print !x}' )

	${software}/run_savi.py --verbose \
	     --input ${outputdir}/${SGE_TASK_ID}.vcf.bgz \
	     --name savi_${SGE_TASK_ID} \
	     --sample ${compsamp} \
	     --prior ${priorstring} \
	     --saviprecision ${savi_precision} \
	     --keepfreqfile ${keepfreqfile} \
	     --outputdir ${outputdir}/savi_${SGE_TASK_ID}/out

	if [ $filter -eq 1 ]; then 

		# FIX THIS FOR NON-COMPARISON
		# vcffilter -f "PD21_L > 0" ${outputdir}/savi_${SGE_TASK_ID}/out/finalsavi.vcf.bgz | sed 's|chrM|chrMT|g' > ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf     
		# vcffilter -f "PD21_L > 0" ${outputdir}/savi_${SGE_TASK_ID}/out/finalsavi.vcf.bgz > ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf     
		# make a string for the vcffilter - e.g., if $compsamp is 2:1,3:1,5 then $filterstring will be PD21_L > 0 | PD31_L > 0
		filterstring=$( echo $compsamp | perl -ne '{chomp($line = $_); @myarr = split(/,/, $line); foreach my $elt (@myarr) {if ($elt =~ s/://) {print "PD".$elt."_L > 0 | "}}}' | sed 's|\| $||' )
		vcffilter -f "${filterstring}" ${outputdir}/savi_${SGE_TASK_ID}/out/finalsavi.vcf.bgz > ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf 

		# index result (too small? not nec?)
		# bgzip ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf
		# mv ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf.gz ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf.bgz
		# tabix -p vcf ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf.bgz
	else
		zcat ${outputdir}/savi_${SGE_TASK_ID}/out/freqsavi.vcf.bgz > ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf     
		rm ${outputdir}/savi_${SGE_TASK_ID}/out/freqsavi.vcf.bgz ${outputdir}/savi_${SGE_TASK_ID}/out/freqsavi.vcf.bgz.tbi
	fi

	# if empty
	if [ $( cat ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf | sed '/^#/d' | head | wc -l ) == 0 ]; then
		echo "[halt script] SAVI vcf file is empty"
		exit 0;
	fi

	echo "[date 4] "`date`;
fi

### SnpEff
outputdir2=${outputdir}/ann_${SGE_TASK_ID}
if [[ $stepstr == *5* ]]; then
	echo "[STEP5] annotation"
	mkdir -p ${outputdir2}

	java -Xmx${java_memory}G -jar $( which snpEff.jar ) -c $( which snpEff.config ) $anngenome -noStats -v -lof  -canon -no-downstream -no-intergenic -no-intron -no-upstream -no-utr ${outputdir}/savi_${SGE_TASK_ID}/out/finalfilter.vcf > ${outputdir2}/tmp_${SGE_TASK_ID}.eff_0.vcf

	# loop through all annotating vcfs
	myindex=1
	for i in $( echo ${annvcf} | sed 's|,| |g' ); do 
		java -Xmx${java_memory}G -jar $( which SnpSift.jar ) annotate ${i} -v ${outputdir2}/tmp_${SGE_TASK_ID}.eff_$(($myindex - 1)).vcf > ${outputdir2}/tmp_${SGE_TASK_ID}.eff_${myindex}.vcf
		if [ $debug -eq 0 ]; then 
			rm ${outputdir2}/tmp_${SGE_TASK_ID}.eff_$(($myindex - 1)).vcf
		fi
		myindex=$(($myindex + 1))
	done

	echo "[date 5] "`date`;
fi

# Copy Number Variation
outputdir3=${outputdir}/cnv_${SGE_TASK_ID}
if [[ $stepstr == *6* ]]; then
	echo "[STEP6] copy number variation"
	### CNV and LOH
	# THIS ONLY WORKS FOR 2 SAMPLES - PAIRED TUMOR NORMAL DATA

	mkdir -p ${outputdir3}/prior

	cat ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt | java -Xmx${java_memory}G \
	    -jar $( which VarScan.v2.3.6.jar ) copynumber \
	    -mpileup \
	    ${outputdir3}/${SGE_TASK_ID}.tumor

	cat ${outputdir3}/${SGE_TASK_ID}.tumor.copynumber | sed '1d' | awk 'BEGIN{OFS="\t"}{print "1","40",int($6+0.5),int($5+$6+0.5)}' > ${outputdir3}/${SGE_TASK_ID}.qvt.txt
	${software}/make_prior.py -n ${SGE_TASK_ID}.tumor -j 10 -p ${software}/savi/unif_prior01 -o ${outputdir3}/prior --qvt ${outputdir3}/${SGE_TASK_ID}.qvt.txt

	mymean=$( cat ${outputdir3}/${SGE_TASK_ID}.qvt.txt | awk 'BEGIN{x=0}{x=x+$3/$4;}END{print x/NR}' )

	cat ${outputdir3}/${SGE_TASK_ID}.qvt.txt | ${software}/savi/savi_poster -p ${outputdir3}/prior/prior | ${software}/savi/savi_comp 2 ${mymean} \
	    | awk 'BEGIN{print "cnv_call\tcnv_call_posterior"}{print}' > ${outputdir3}/${SGE_TASK_ID}.cnv.tmp.txt

	paste ${outputdir3}/${SGE_TASK_ID}.tumor.copynumber ${outputdir3}/${SGE_TASK_ID}.cnv.tmp.txt > ${outputdir3}/${SGE_TASK_ID}.cnv.txt

	# Loss of Heterozygosity
	zcat $outputdir/${SGE_TASK_ID}.vcf.bgz | ${software}/savi/make_qvt -1 -2s 2,1 \
	   | ${software}/savi/savi_poster -pd ${outputdir}/savi_${SGE_TASK_ID}/prior2/prior ${outputdir}/savi_${SGE_TASK_ID}/prior1/prior \
	   | ${software}/savi/savi_comp 2 0 \
	   | awk 'BEGIN{print "loh_call\tloh_call_posterior"}{print}' > ${outputdir3}/${SGE_TASK_ID}.loh.txt

	echo "[date 6] "`date`;
fi

echo "[STEP CLEANUP] clean up"
# if not debug
if [ $debug -eq 0 ]; then 
	rm -f ${outputdir}/tmp_mpile.${SGE_TASK_ID}.txt 
	rm -f ${outputdir3}/${SGE_TASK_ID}.qvt.txt
	rm -f ${outputdir3}/${SGE_TASK_ID}.cnv.tmp.txt
	rm -f ${outputdir3}/${SGE_TASK_ID}.tumor.copynumber
fi

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
echo [End] `date`
