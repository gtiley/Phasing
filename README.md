----Using the script----
Try this to see options:
perl phaseLoci.pl

Important note - you will run the script twice:
	  a) First in --runmode 0 or 1
	  b) Second in --runmode 2

Important note - you will need a copy of the H-PoPGv0.2.0.jar file
https://github.com/MinzhuXie/H-PoPG

Important note - this is under active development, it works but we have several things planned in the near future
	  1) More user-friendly
	  2) GATK filter options
	  3) Better summary statistics
	  4) Ploidy estimation 

----Some notes on configuring data and folders----
1: Naming of Fastq files
   Fastq files follow the following naming rules:
   	 a) Only paired-end data allowed
	 b) Reads should be named as
   	    <Individual ID>.R1.<Fastq File Extension>
	    <Individual ID>.R2.<Fastq File Extention>
	    where <Individual ID> = The individual name specified in the ploidy file
	    and <Fastq File Extension> = Whatever you want; It does not matter if named *.fq, *.fastq, *.fq.gz, etc
   
   Fastq files are assumed to be pre-processed for quality and adapter removal. We do not integrate such tools here as some would argue that letting the soft-clipping in BWA is a better approach and GATK deals with quality explicitly.

2: The Ploidy File
   This is where Individual IDs and their ploidy levels are specified. It has the following rules:
   	a)At least two columns are present
	b)The first column is Individual ID
	c)The second column is the ploidy level represented by an integer
	d)More columns are allowed with any metadata you would like for your own purposes
	e)Columns are seperated by whitespaces, do not format as comma-seperated

3: The Reference Fastas
   Reference fasta files have the following rules:
   	a) There is a sinlge fasta file per locus that contains reference sequences for all individuals
	b) The reference fasta is named <Locus Name>.fasta
	   i) No spaces allowed in locus name
	   ii) The locus name here will be the locus name of the output fasta with phased sequences
	c) The fasta files are not interleaved and assume no line breaks in the sequence data
	d) The fasta headers are assumed to match the Individual ID (i.e. ><Individual ID>)
	e) Other information can follow the Individual ID in the fasta header, but will be ignored      

4: The Template File
   There is a file called template.sh to help distribute jobs on a cluster
        a) Make the necessary changes to the scheuduler directives for your account
        b) I was on a cluster with SLURM when making this - you will need to edit for PBS or SGE accordingly
        c) If you run the script in runmode 1, all of the commands are pasted below what you have in the template.sh file - all you need to do is configure the directives and paths to software here if they are not specified in the control file
	d) The template file must always be provided as an argument, even if using runmode 0 or 2
	e) The template file can be renamed if you like