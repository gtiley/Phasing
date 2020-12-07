# Using the script
## Try this to see options:
```perl phaseLoci.pl```

## Important note - you will run the script twice
	  1. First in --runmode 0 or 1
	  	* 0 will run the phasing pipeline one each individual, one at a time, on a single processor. This might be helpful when analyzing a small number of individuals on a local machine. Otherwise, I recommend 1 so each individual can be distributed on a single processesor when on a cluster. All of your necessary cluster configuration needs to happen in *template.sh*.
		* 1 will distribute 
	  2. Second in --runmode 2
	  	* 2 generates the per-locus fasta output and summary statistics. This step is fast and happens on a single processor. Here the genotype <0 || 1> option comes into effect.
		* The genotype option affects how variants with ambiguous phases are handeled. When multiple haplotype blocks are recovered for a locus, we retain the phasing of the block with the most variants only. 0 causes the others to be replaced with "N" while 1 causes these unphased variants to be replaced with there IUPAC codes. There may be analyses where one option is more favorable than the other, so we make both possibilities available here.  

### Important note - you will need a copy of the H-PoPGv0.2.0.jar file from https://github.com/MinzhuXie/H-PoPG

## Important note - this is under active development, it works but we have several things planned in the near future
	  1. Data preprocessing
	  2. GATK filter options
	  3. Ploidy estimation

# Some notes on configuring data and folders
## Naming of Fastq files
Fastq files follow the following naming rules:
* Only paired-end data allowed
* Reads should be named as
	+ &lt;Individual ID&gt;.R1.&lt;Fastq File Extension&gt;
	+ &lt;Individual ID&gt;.R2.&lt;Fastq File Extention&gt;
	+ where &lt;Individual ID&gt; = The individual name specified in the ploidy file
	+ and &lt;Fastq File Extension&gt; = Whatever you want; It does not matter if named .fq, .fastq, .fq.gz, etc
Fastq files are assumed to be pre-processed for quality and adapter removal. We do not integrate such tools here as some would argue that the soft-clipping in BWA is a better approach and GATK deals with quality explicitly.

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
