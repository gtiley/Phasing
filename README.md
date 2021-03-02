# Using the script
## Try this to see options:
```perl phaseLoci.pl```

## There are four input options and all are required
```perl phaseLoci.pl --controlFile phase.ctl --runmode 1 --template template.sh --genotype 0```
1. The control file (phase.ctl) is used to configure the paths to other software, your input data, and where you want your output.
2. The runmode flag is used to determine whether you are running the pipeline serially on perhaps your personal computer or wanting to distribute jobs on a cluster. See below for more details.
3. The template file (template.sh) has the basic directives you will use for running on a cluster.
4. The genotype flag determines how phasing ambiguity is handeled. See below for details.

## Important Note - you will run the script twice
1. First in --runmode 0 or 1
+ 0 will run the phasing pipeline one each individual, one at a time, on a single processor. This might be helpful when analyzing a small number of individuals on a local machine. Otherwise, I recommend 1 so each individual can be distributed on a single processesor when on a cluster. All of your necessary cluster configuration needs to happen in *template.sh*.
+ 1 will distribute 
2. Second in --runmode 2
+ 2 generates the per-locus fasta output and summary statistics. This step is fast and happens on a single processor. Here the genotype <0 || 1> option comes into effect.
+ The genotype option affects how variants with ambiguous phases are handeled. When multiple haplotype blocks are recovered for a locus, we retain the phasing of the block with the most variants only. 0 causes the others to be replaced with "N" while 1 causes these unphased variants to be replaced with there IUPAC codes. There may be analyses where one option is more favorable than the other, so we make both possibilities available here.  

## Important Note - you will to install a few software on your computer or cluster
Please cite and credit the authors of all of the important bits that are glued together here.
* BWA (Li et al. 2009a)
* The samtools/htslib library (Li et al. 2009b)
* GATK (McKenna at al. 2012)
* HPoPG (Xie et al. 2016)

## Important Note - this is under active development, it works but we have several things planned in the near future
1. Data preprocessing
2. GATK filter options
+ This is maybe the most immediate need. Filters might require fine-tuning depending on the study. Currently the filtering expressions are hardcoded as described in the manuscript. If you want edit filters, this should be intuitive when inspecting line 275 for runmode=0 or line 371 for runmode=1.
3. Ploidy estimation
4. Incorporating joint genotyping with a single reference

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

## The Template File
There is a file called template.sh to help distribute jobs on a cluster
* Make the necessary changes to the scheuduler directives for your account
* I was on a cluster with SLURM when making this - you will need to edit for PBS or SGE accordingly
* If you run the script in runmode 1, all of the commands are pasted below what you have in the template.sh file - all you need to do is configure the directives and paths to software here if they are not specified in the control file
* The template file must always be provided as an argument, even if using runmode 0 or 2
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
Second, there has apparently been a flurry of phasing algorithms developed for polyploids in the recent years after my colleagues and I began working on this pipeline and our ideas. I will try to investigate and compare some of them (Moeinzadeh et al. 2020; ) in the future 


# References
* Li H, Durbin R. 2009a. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754-1760.
* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. 2009b. The sequence alignment/map format and SAMtools. Bioinformatics 25:2078-2079. 
* McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, et al. 2010. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20:1297-1303.
* Xie M, Wu Q, Wang J, Jiang T. 2016. H-PoP and H-PoPG: heuristic partitioning algorithms for single individual haplotyping of polyploids. Bioinformatics 32:3735-3744.
* Gerard D, Ferrão LFV, Garcia AAF, Stephens M. 2018. Genotyping polyploids from messy sequencing data. Genetics 210:789-807.
* Moeinzadeh M-H, Yang J, Muzychenko E, Gallone G, Heller D, Reinert K, Haas S, Vingron M. 2020. Ranbow: A fast and accurate method for polyploid haplotype reconstruction. PLoS Comp Biol. 16:e1007843.
* Schrinner SD, Mari RS, Ebler J, Rautiainen M, Seillier L, Reimer JJ, Usadel B, Marschall T. 2020. Haplotype threading: an accurate polyploid phasing from long reads. Genome Biol. 21:252.
