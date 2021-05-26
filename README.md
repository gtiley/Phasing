# Using the script
## Try this to see options:
```perl PATE.pl```

## There are four input options and all are required
```perl PATE.pl --controlFile PATE.ctl --runmode cluster --template template.sh --genotype consensus```
1. The control file (PATE.ctl) is used to configure the paths to other software, your input data, and where you want your output.
2. The runmode flag is used to determine whether you are running the pipeline serially on perhaps your personal computer or wanting to distribute jobs on a cluster. See below for more details.
3. The template file (template.sh) has the basic directives you will use for running on a cluster.
4. The genotype flag determines how phasing ambiguity is handeled. See below for details.

## Important Note - you will run the script twice
1. First in --runmode ```serial``` or ```cluster```
+ ```serial``` will run PATÉ on each individual, one at a time, on a single processor. This might be helpful when analyzing a small number of individuals on a local machine. Otherwise, I recommend 1 so each individual can be distributed on a single processesor when on a cluster. All of your necessary cluster configuration needs to happen in *template.sh*.
+ ```cluster``` will distribute 
2. Second in --runmode ```alleles```
+ ```alleles``` generates the per-locus fasta output and summary statistics. This step is fast and happens on a single processor. Here the genotype < ```consensus``` || ```iupac``` > option comes into effect.
+ The genotype option affects how variants with ambiguous phases are handeled. When multiple haplotype blocks are recovered for a locus, we retain the phasing of the block with the most variants only. ```consensus``` causes the others to be replaced with "N" while ```iupac``` causes these unphased variants to be replaced with there IUPAC codes. There may be analyses where one option is more favorable than the other, so we make both possibilities available here.  

## Important Note - you will to install a few software on your computer or cluster
Please cite and credit the authors of all of the important bits that are glued together here.
* BWA (Li et al. 2009a)
* The samtools/htslib library (Li et al. 2009b)
* GATK (McKenna at al. 2012)
* HPoPG (Xie et al. 2016)

## Important Note - this is under active development, it works but we have several things planned in the near future
1. Ploidy estimation
2. Parent assignment
3. Determining allo vs. autopolyploidy
4. Incorporating joint genotyping with a single reference

## Explanation of the control file options
### There are several input and output folders and files to keep track of
* REF = input folder reference fasta files
* GENOTYPE\_OUT = output folder for genotyping files
* PHASE\_OUT = output folder for phased sequences, but split by individual
* FASTA\_OUT = output folder for fasta files you want to use
+ PHASED - these are the fasta files with all phased alleles per-locus
+ GENOTYPE - these fasta files have the unphased genotype sequences
* IUPAC\_OUT = output folder for phased sequences, but split by individual
* SUMMARYSTATS\_OUT = output folder with all of the phasing summary statisitcs
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
Second, there has apparently been a flurry of phasing algorithms developed for polyploids in the recent years after my colleagues and I began working on this pipeline and our ideas. I will try to investigate and compare some of them (e.g. Moeinzadeh et al. 2020) in the future 


# References
* Li H, Durbin R. 2009a. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754-1760.
* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. 2009b. The sequence alignment/map format and SAMtools. Bioinformatics 25:2078-2079. 
* McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, et al. 2010. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20:1297-1303.
* Xie M, Wu Q, Wang J, Jiang T. 2016. H-PoP and H-PoPG: heuristic partitioning algorithms for single individual haplotyping of polyploids. Bioinformatics 32:3735-3744.
* Gerard D, Ferrão LFV, Garcia AAF, Stephens M. 2018. Genotyping polyploids from messy sequencing data. Genetics 210:789-807.
* Moeinzadeh M-H, Yang J, Muzychenko E, Gallone G, Heller D, Reinert K, Haas S, Vingron M. 2020. Ranbow: A fast and accurate method for polyploid haplotype reconstruction. PLoS Comp Biol. 16:e1007843.
* Schrinner SD, Mari RS, Ebler J, Rautiainen M, Seillier L, Reimer JJ, Usadel B, Marschall T. 2020. Haplotype threading: an accurate polyploid phasing from long reads. Genome Biol. 21:252.
