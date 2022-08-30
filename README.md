# Using the script
## Try this to see options:
```perl PATE.pl```

## There are four input options and all are required
```perl PATE.pl --controlFile PATE.ctl --runmode species --template template.sh --genotype consensus```
1. The control file (PATE.ctl) is used to configure the paths to other software, your input data, and where you want your output.
2. The runmode flag is used to determine if genotyping will be done jointly or not. There may be different conditions where one is desirable over the other. See below for more details.
3. The template file (template.sh) has the basic directives you will use for running on a cluster. Or if running on your local computer, you can set environment variables at run time.
4. The genotype flag determines how phasing ambiguity is handeled. See below for details.

## Important Note - you will run the script twice
1. First in --runmode ```species``` or ```population1```
	+ ```species``` will run GATK *HaplotypeCaller* on each individual, that is realigning short reads to the consensus assemblies to get genotypes based on ploidy levels specified in the *ploidy file*.
	+ ```population1``` will run GATK on each individual in gVCF mode and set up the script for joint genotyping based on the specified ploidy levels. You will have to submit the joint genotyping job separately after the gVCF files are generated. The reference individual is specified in the *control file*.
2. Second in --runmode ```alleles``` or ```population2```
	+ ```alleles``` generates the per-locus fasta output and summary statistics based on the individual genotype data. This step is fast and happens on a single processor. Here the genotype < ```consensus``` || ```iupac``` > option comes into effect.
	+```population2``` is similar to *alleles* but will split the multisample VCF by individual to pass to *H-PoPG* to handle the phasing. The same genotype flags apply.
	+ The genotype option affects how variants with ambiguous phases are handeled. When multiple haplotype blocks are recovered for a locus, we retain the phasing of the block with the most variants only. ```consensus``` causes the others to be replaced with "N" while ```iupac``` causes these unphased variants to be replaced with there IUPAC codes. There may be analyses where one option is more favorable than the other, so we make both possibilities available here.  

## Important Note - you will to install a few software on your computer or cluster
Please cite and credit the authors of all of the important bits that are glued together here.
* BWA (Li et al. 2009a)
* The samtools/htslib library (Li et al. 2009b)
* GATK (McKenna at al. 2012)
* HPoPG (Xie et al. 2016)

## Experimental feature
It is possible to estimate ploidy directly for the distributions of ratios with reads that carry the reference or alternate allele (Weiß et al. 2018; Viruel 2019). This is achieved by the --runmode ```estPloidy``` option and specifying an appropriate reference in the control file. All individuals are genotyped as a diploid and mixture models (Tiley et al. 2018) are used to determine the ploidy. I am skeptical that enough information is available in target enrichment data to achieve this task and an outgroup that can successfully polarize the alternate alleles are needed. Nevertheless, the option exists and may have value for some cases.

## Important Note - this is under active development, it works but we have several things planned in the near future
1. Imputation
2. SFS estimation
2. Parent assignment
3. Determining allo vs. autopolyploidy

## Explanation of the control file options
### There are several input and output folders and files to keep track of
* PHASING_ROOT = the root directory that you clone or download from github
* REF = input folder reference fasta files
* GENOTYPE\_OUT = output folder for genotyping files
* PHASE\_OUT = output folder for phased sequences, but split by individual
* FASTA\_OUT = output folder for fasta files you want to use
	+ PHASED - these are the fasta files with all phased alleles per-locus
	+ GENOTYPE - these fasta files have the unphased genotype sequences
	+ PICKONE - these are the fasta files with only one haplotypic sequence available per individual per locus. These data may be preferred for phylogenetic network analyses.
* IUPAC\_OUT = output folder for phased sequences, but split by individual
* SUMMARYSTATS\_OUT = output folder with all of the phasing summary statistics
* ESTPLOIDY\_OUT = output folder with all genotyping and mixture model results along with a new ploidy file
* FQ = input folder of fastq files
* PLOIDY = input ploidy file (see below)

### The software dependencies
* BWA = path to bwa
* PICARD = path to picard
* GATK = path to gatk
* SAMTOOLS = path to samtools
* BAMTOOLS = path to bamtools
* BGZIP = path to bgzip from htslib 
* TABIX = path to tabix from htslib
* HPOPG = path to H-PoPG jar file
* SCHEDULER = The submission command for your scheduler (SLURM uses sbatch)

### We make the VCF filtering options for GATK an option that the user can change. Each line is a separate filter and has the specific formatting of "<FILTER>" <whitespace> "<FILTER NAME>"
### The following are some reasonable defaults, but will note always be optimal. Consult VCF filtering annotations before altering.
* GATK\_FILTER\_EXPRESSION = "QD < 2.0" "QD\_lt2"
* GATK\_FILTER\_EXPRESSION = "FS > 60.0" "FS\_gt60"
* GATK\_FILTER\_EXPRESSION = "MQ < 40.0" "MQ\_lt40"
* GATK\_FILTER\_EXPRESSION = "ReadPosRankSum < -8.0" "ReadPosRankSum\_ltm8"
* GATK\_FILTER\_EXPRESSION = "AF < 0.05" "AF\_lt05"
* GATK\_FILTER\_EXPRESSION = "AF > 0.95" "AF\_gt95"
* GATK\_FILTER\_EXPRESSION = "DP < 10" "DP\_lt10"

### Some additional options
* UNIQUE\_ONLY = <0 || 1> determines if only unique haplotypes are output following phasing by H-PoPG (0=no or 1=yes). It is very possible to have multiple identical haplotypes per locus.
* REMOVE\_READ\_DUPLICATES = <0 || 1> = should Picard be used to remove pcr duplicates. Leave at 0=no for all target enrichment studies but this can be turned on for non-enriched libraries.
* OUTPUT\_EXPECTED\_DOSAGE = <0 || 1> = should the expected number of sequences per individual be output even in the absence of variants. This might be needed when calculating allele frequencies.

To generate the output comparable to Tiley et al. (2021) or Crowl et al. (2022), use the following options:
* UNIQUE\_ONLY = 1
* REMOVE\_READ\_DUPLICATES = 0
* OUTPUT\_EXPECTED\_DOSAGE = 0

# Some notes on configuring data and folders
## Naming of Fastq Files
Fastq files follow the following naming rules:
* Only paired-end data allowed
* Reads should be named as
	+ &lt;Individual ID&gt;.R1.&lt;Fastq File Extension&gt;
	+ &lt;Individual ID&gt;.R2.&lt;Fastq File Extention&gt;
	+ where &lt;Individual ID&gt; = The individual name specified in the ploidy file
	+ and &lt;Fastq File Extension&gt; = Whatever you want; It does not matter if named .fq, .fastq, .fq.gz, etc
* Fastq files are assumed to be pre-processed for quality and adapter removal. We do not integrate such tools here as some would argue that the soft-clipping in BWA is a better approach and GATK deals with quality explicitly.
* A helper script is available to format the fastq names for you, please see ```helperScripts/PATE_formatInput.pl```

## The Ploidy File
This is where Individual IDs and their ploidy levels are specified. It has the following rules:
* At least two columns are present
* The first column is Individual ID
* The second column is the ploidy level represented by an integer
* More columns are allowed with any metadata you would like for your own purposes
* Columns are seperated by whitespaces, do not format as comma-seperated

## The Reference Fastas
Reference fasta files have the following rules:
* There is a sinlge fasta file per locus that contains reference sequences for all individuals
* The reference fasta is named &lt;Locus Name&gt;.fasta
	+ No spaces allowed in locus name
	+ The locus name here will be the locus name of the output fasta with phased sequences
* The fasta files are not interleaved and assume no line breaks in the sequence data
* The fasta headers are assumed to match the Individual ID (i.e. &gt;&lt;Individual ID&gt;)
* Other information can follow the Individual ID in the fasta header, but will be ignored
* A helper script is available to format the fasta files for you, please see ```helperScripts/PATE_formatInput.pl```      

## The Template File
There is a file called template.sh to help distribute jobs on a cluster
* Make the necessary changes to the scheuduler directives for your account
* I was on a cluster with SLURM when making this - you will need to edit for PBS or SGE accordingly
* If you run the script in runmode ```cluster```, all of the commands are pasted below what you have in the template.sh file - all you need to do is configure the directives and paths to software here if they are not specified in the control file
* The template file must always be provided as an argument, even if using runmode ```serial``` or ```alleles```
* The template file can be renamed if you like

# Explanation of Summary Statistics
An output after runmode=2 is ```averagePhasingStats.txt```, which contains a few numbers aggregated over all loci. Here is a description of each column.
1. INDIVIDUAL - individual ID
2. NLOCI - number of loci assembled by hybpiper
3. NVARLOCI - number of loci with at least 1 variant        
4. NINVLOCI - number of loci with no variants        
5. NPHASELOCI - number of loci phased (it is possible to have unphased loci with variants, in the case when there is 1 or more blocks with only 1 variant)      
6. AVG_LENGTH - average locus length      
7. AVG_NVAR - average number of variants (includes invariable loci - averaged across the total number of loci)        
8. AVG_HET - average heterozygosity (per-base heterozygosity averaged over all loci - including the invariable ones) 
9. AVG_NBLOCKS - average number of blocks per phased locus (thus, excluding loci without phasing information. Many of the Dryopteris individuals were close to 1 or even 1 in some cases. I think this is incredibly useful for evaluating how good the phasing is working.)
10. AVG_LONGESTBL - Average number of variants in the longest phasing block. (thus not including loci without phasing information) 
	+ 1 - number of loci that were completely homozygous       
	+ 2 - number of loci with two phased alleles	      
	+ 3       .
	+ 4       .
	+ 5       .
	+ 6 - number of loci with six phased alleles (six was the max in Dryopteris. If you have higher ploidies, this will automatically go higher because it is based off of the ploidy file. If you have octoploids for example, this will then go to 8)

# Opinions
There are several technical issues compounded in the existing pipeline and I view this as a starting point for enabling some interesting analyses of polyploid complexes. 
First, the genotyping problem in polyploids has a lot of uncertainty and I recommend reading Gerard et al. (2018) to better appreciate the problems. Comparisons of genotypers are needed in the future, but we do our best to filter out errors from the final set of variants.
Second, there has apparently been a flurry of phasing algorithms developed for polyploids in the recent years after my colleagues and I began working on this pipeline and our ideas. I will try to investigate and compare some of them (e.g. Moeinzadeh et al. 2020) in the future.


# References
* Crowl AA, Fritsch PW, Tiley GP, Lynch NP, Ranney TG, Ashrafi H, Manos PS. 2022. A First Complete Phylogenomic Hypothesis for Diploid Blueberries (Vaccinium section Cyanococcus). American Journal of Botany. In press.
* Li H, Durbin R. 2009a. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754-1760.
* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. 2009b. The sequence alignment/map format and SAMtools. Bioinformatics 25:2078-2079. 
* McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, et al. 2010. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20:1297-1303.
* Xie M, Wu Q, Wang J, Jiang T. 2016. H-PoP and H-PoPG: heuristic partitioning algorithms for single individual haplotyping of polyploids. Bioinformatics 32:3735-3744.
* Gerard D, Ferrão LFV, Garcia AAF, Stephens M. 2018. Genotyping polyploids from messy sequencing data. Genetics 210:789-807.
* Moeinzadeh M-H, Yang J, Muzychenko E, Gallone G, Heller D, Reinert K, Haas S, Vingron M. 2020. Ranbow: A fast and accurate method for polyploid haplotype reconstruction. PLoS Comp Biol. 16:e1007843.
* Schrinner SD, Mari RS, Ebler J, Rautiainen M, Seillier L, Reimer JJ, Usadel B, Marschall T. 2020. Haplotype threading: an accurate polyploid phasing from long reads. Genome Biol. 21:252.
* Tiley GP, Barker MS, Burleigh JG. 2018. Assessing the Performance of *Ks* Plots for Detecting Ancient Whole Genome Duplications. Genome Biology and Evolution 10:2882-2898.
* Tiley GP, Crowl AA, Manos PS, Sessa EB, Solís-Lemus C, Yoder AD, Burleigh JG. 2021. Phasing alleles improves network inference with allopolyploids. bioRxiv doi: https://doi.org/10.1101/2021.05.04.442457
* Viruel J, Conejero M, Hidalgo O, Pokorny L, Powell RF, Forest F, Kantar MB, Soto Gomez M, Graham SW, Gravendeel B, Wilkin P, Leitch IJ. 2019. A Target Capture-Based Method to Estimate Ploidy from Herbarium Specimens. Front. Plant Sci. 10:937.
* Weiß CL, Pais M, Cano LM, Kamoun S, Burbano HA. 2018. nQuire: a statistical framework for ploidy estimation using next generation sequencing. BMC Bioinformatics 19:122.
