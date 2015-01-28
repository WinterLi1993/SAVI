## SAVI
SAVI - *statistical algorithm for variant identification* - is a program for calling variants in sequencing data written by Vladimir Trifonov with contributions from Oliver Elliott, Jiguang Wang, Alex Penson, and Hossein Khiabanian.

**Introduction**
Why run SAVI? In a word, SAVI is for finding needles in haystacks. It can boil a very large dataset into a small list of mutations. In bioinformatics, *calling variants* can mean two different things---(1) simply enumerating all differences between some sequence data and a reference; and (2) determining which of those differences are significant and not likely to be error. SAVI does the later, while a program like Samtools mpileup will do the former. In practice, SAVI is a way of sifting through large amounts of data to pull out the significant mutations using a Bayesian probability model. A common use case is identifying deleterious mutations in cancer, given normal and tumor sequence data---an often-encountered problem in bioinformatics. The output of SAVI is a list of candidate genomic alterations each annotated with a probability of how likely it is to be real. 

SAVI works with standard bioinformatic file formats. As input, the SAVI pipeline takes a bam files and it produces a vcf file as output. In the output vcf file, the SAVI probabilities are added in the INFO field. 

If you're interested in the mathematical underpinings of SAVI, you can read about it in [this BMC Systems Biology paper](https://www.biomedcentral.com/1752-0509/7/S2/S2).

**Dependencies**

The following programs must be in your PATH:

- Samtools
- SnpEff (i.e., which snpEff.jar must return a path)
- [tabix](http://samtools.sourceforge.net/tabix.shtml)
- [bgzip](http://samtools.sourceforge.net/tabix.shtml)
- [vcflib](https://github.com/ekg/vcflib)
- (VarScan, if you want to run the Copy Number Variation and Loss of Heterozygosity Computation step)

**Installation**

Installing SAVI is simple. First, clone this repo. Then `make` the binaries in the SAVI/savi directory:

```
cd SAVI/savi
make
```

**Additional Files**

SAVI needs a reference fasta file (whatever you mapped your bams to) and its faidx index.
It is also helpful to use various vcf files to provide additional annotation.
You can download all of this on the Rabadan lab homepage [here](http://rabadan.c2b2.columbia.edu/public/savi_resources/).
You'll find the following files. A human reference, hg19:

- hg19_chr.fold.25.fa
- hg19_chr.fold.25.fa.fai

And various annotating vcfs:

- dbSnp138.vcf - dbSnp 138
- cbio.fix.sort.vcf - cBio variants
- CosmicVariants_v66_20130725.vcf - Cosmic variants
- 219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf - Rabadan Lab supernormal
- meganormal186TCGA.fix.sort.vcf - Rabadan Lab TCGA supernormal

**Workflow**

The SAVI pipeline is organized into 6 main steps. The steps are:

1. Format conversion: bams to mpileup (input: bam files; output: pileup file)
2. Format conversion: mpileup to multiallelic vcf (input: pileup file; output: vcf file)
3. Savi: Make Prior (input: vcf file; output: prior files)
4. Savi: Run Savi (input: vcf file + prior files; output: filtered vcf file)
5. SnpEff Annotation (input: filtered vcf file; output: filtered, annotated vcf file)
6. Copy Number Variation and Loss of Heterozygosity Computation (optional for the Tumor-Normal paired samples only) (input: pileup; output: text file)

**Usage Examples**

The most common use case for SAVI is paired normal tumor bam files (mapped to the same reference, of course) where the goal is to find the significant mutations in the tumor sample.
Before we run SAVI, we need to make sure of the following:

- the reference fasta we're providing to SAVI is the same one to which we've mapped our bam files
- we've indexed the reference fasta with faidx
- we've sorted our bams files via `samtools sort`
- we've indexed our sorted bam files via `samtools index`

Got it? Good! Here's an example command running SAVI **for chromosome 1**:

```
SAVI/run_pipeline.sh --bam normal.bam,tumor.bam --ref savi_resources/hg19_chr.fold.25.fa --memory 4 --compsamp 2:1 --outputdir outputdir/samplename/chr1 --region chr1 --ann-vcf savi_resources/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf,savi_resources/cbio.fix.sort.vcf,savi_resources/CosmicVariants_v66_20130725.vcf,savi_resources/dbSnp138.vcf,savi_resources/meganormal186TCGA.fix.sort.vcf
```

In practice, we'd want to run SAVI for every chromosome in a loop:

```
for i in {1..22} X Y M; do
	# command here
done
```

Another use case is a tumor-only bam file.
Here's a sample command again **for chromosome 1**:

```
SAVI/run_pipeline.sh --bam tumor.bam --ref savi_resources/hg19_chr.fold.25.fa --memory 4 --compsamp 1 --nofilter --outputdir outputdir/samplename/chr1 --region chr1 --ann-vcf savi_resources/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf,savi_resources/cbio.fix.sort.vcf,savi_resources/CosmicVariants_v66_20130725.vcf,savi_resources/dbSnp138.vcf,savi_resources/meganormal186TCGA.fix.sort.vcf
```

In this case, SAVI won't be able to call somatic variants.
Instead, it will merely tell us what it thinks are present---a much longer list than if we had the normal sample for comparison.

**Common Problems**

- Are your bams sorted by position? - samtools mpileup requires "position sorted alignment files."
- Did you index your bams to produce bai files? - Samtools won't be able to extract regions of your bam files if they have not been indexed.
- Are both the tumor and normal bam files mapped to the same reference?

**More Detailed Documentation**

For more detailed documentation, refer to [this page](goo.gl/LM9iMN).
