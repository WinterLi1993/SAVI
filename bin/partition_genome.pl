#!/usr/bin/env perl

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  About:
#
#  Partition a File of Genomic Ranges
#
#  Useage example:
#  cat genome_range.txt | $0 1000 
#  
#  Args:
#  first argument - partition size 
#
#  Example
#
#  \$ cat hg19.chrom_24chrs.sizes | head -2
#  chr1    249250621
#  chr2    243199373
#  \$ cat hg19.chrom_24chrs.sizes | head -2 | ./partition_genome.pl 50000000 
#  chr1:1-50000000
#  chr1:50000001-100000000
#  chr1:100000001-150000000
#  chr1:150000001-200000000
#  chr1:200000001-249250621
#  chr2:1-50000000
#  chr2:50000001-100000000
#  chr2:100000001-150000000
#  chr2:150000001-200000000
#  chr2:200000001-243199373
#
####################################################################################################################

_EOUSAGE_

#my $arg1 = shift;
#my $arg2 = shift;

my $arg1 = $ARGV[0];

if ($arg1 eq "-h" or $arg1 eq "--h" or $arg1 eq "-help" or $arg1 eq "--help" or scalar(@ARGV) == 0)
{
	print $usage;
	exit;
}


while (<STDIN>)
{	
	chomp($_); 
	$myline = $_; 
	@a=split(/\s/, $_); 
	$mychr=$a[0]; 
	# $mydiv = 50000000; 
	$mydiv = $arg1;
	$mynum = $a[1]; 
	$myquo = int($mynum/$mydiv); 
	$myremainder = $mynum%$mydiv; 
	if (1==0) {print $myline."\t".$myquo."\t".$myremainder."\n"}; 
	$range2=0; 
	for (my $i = 1; $i <= $myquo; $i++) 
	{
		$range1=$i*$mydiv-$mydiv+1; $range2=$i*$mydiv; print $mychr.":".$range1."-".$range2."\n";
	} 
	if ($myremainder>0) 
	{
		$newrange1=$range2 + 1; 
		$newrange2=$range2+$myremainder; 
		print $mychr.":".$newrange1."-".$newrange2."\n"
	}
}
