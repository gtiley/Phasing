#!/usr/bin/perl -w
use strict;

#----------------------------------------------------------------------------------------#
#30 August 2022
#contact regarding code details: George P. Tiley g.tiley@kew.org
#contact regarding empirical performance: Andy Crowl andrew.crowl@duke.edu
#Genotype target-enrichment loci of known ploidy
#Estimate ploidy from allele balance if unknown
#Phase haplotype sequences for each locus
#----------------------------------------------------------------------------------------#

# Accepted AND necessary commands
my @checkArgs = ("controlFile","runmode","template","genotype");
my %passedArgs = ();
if (scalar(@ARGV) == 0)
{
die "/*--------INPUT PARAMETERS--------*/\n
--controlFile STRING <control file>
--runmode STRING <species or alleles or estPloidy or population1 or population2>
--template STRING <template file for distributing on a cluster>
--genotype STRING <consensus or iupac>

\n/*--------EXAMPLE COMMAND--------*/\n
    perl PATE.pl --controlFile PATE.ctl --runmode species --template template.sh --genotype consensus\n

\n/*--------FLAG OPTIONS--------*/\n
--runmode
	species = genotype each individual against its own de novo assemblies aas a reference with directives given in template.sh
	alleles = use phasing information to make orthologous fasta files with haplotype sequences AFTER --runmode species
	estPloidy = estimate the ploidy for each individual in the ploidy.txt file based on allele balance information when genotying them as diploids. This requires a specified reference.
	population1 = use joint genotyping of all individuals defined in the ploidy.txt file against a single reference at their specified ploidy
	population2 = phase loci for each individual based on output from joint genotyping. The population mode is broken into 2 steps as the amount of time needed for the joint genotyping is highly variable.
	
--genotype
	consensus = Use the original reference sequence for each sample for each locus when there are no variants to phase
	iupac = Use genotyped sequence with IUPAC codes when there are variants, but cannot be phased
	
\n/*--------NOTES--------*/\n
This scripts automates a lot of the genotyping and phasing process, but has to be ran twice.
The first time use --runmode species OR --runmode population. The second time you will use --runmode alleles to create the phased fasta files. If all individuals can be compared to single reference individual, it is possible to estimate the ploidys directly from the sequence data with --runmode estPloidy.
All steps assume a single processor is used or allocated. This takes about 30 minutes for a sample of >100x but <1000x coverage. If you want to change these options you may have to alter the appropriate lines relevant to BWA and GATK.
We have yet to explore effects of different hard filtering strategies through GATK, this might require fine-tuning for low (e.g. <=15x) or very high (>1000x) depth data. Some reasonable defaults (in my opinion GPT) are set in the control file.

\n/*--------Important note affecting output--------*/\n
The control file variable UNIQUE_ONLY is by default set to 0. This creates a scenario where the expected number of alleles are output for every locus for every individual. Thus, there could be identical haplotype sequences output. This can be important for some population genetic analyses. It will cause the number of alleles recovered in the stats files to always match the ploidy level though.\n
Because it can be difficult to know the difference between true homozygosity and false negatives from low read coverage ect, it is possible to only output the unique haplotype sequences. This might be more appropriate for various phylogenetic applications. To only output the unique haplotypes, simply set UNIQUE_ONLY = 1. We used it in earlier iterations of PATÉ but the needs changed with research interests. It did not seem fair to force a user into one scenario versus another, so it is optional.\n\n";
}

elsif (scalar(@ARGV) > 0)
{
    for my $i (0..(scalar(@ARGV) - 1))
    {
		if ($ARGV[$i] eq "--controlFile")
		{
	    	$passedArgs{controlFile} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "--runmode")
		{
	    	$passedArgs{runmode} = $ARGV[$i+1];
	    	if (($passedArgs{runmode} ne "species") and ($passedArgs{runmode} ne "population1") and ($passedArgs{runmode} ne "population2") and ($passedArgs{runmode} ne "alleles") and ($passedArgs{runmode} ne "estPloidy"))
	    	{
	    		die ("ERROR\n\n$passedArgs{runmode} is not a valid runmode option\n\nstopping to fix error\n")
	    	}
		}
		if ($ARGV[$i] eq "--template")
		{
        	$passedArgs{template} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "--genotype")
        {
		    $passedArgs{genotype} = $ARGV[$i+1];
		    if (($passedArgs{genotype} ne "consensus") and ($passedArgs{genotype} ne "iupac"))
	    	{
	    		die ("ERROR\n\n$passedArgs{runmode} is not a valid genotype option\n\nstopping to fix error\n")
	    	}
        }
	}
	foreach my $arg (@checkArgs)
	{
		if (! exists $passedArgs{$arg})
		{
	    	die "/*--------MISSING PARAMETER--------*/\nMissing command line argument: $arg\n\nIf running serially, try --template 0\n\n";
		}
	}
}

my %controlArgs = ();
my %variantFilters = (); 
open FH1,'<',"$passedArgs{controlFile}";
while (<FH1>)
{
	if (/PHASING_ROOT\s+\=\s+(\S+)/)
    {
		$controlArgs{PHASING_ROOT} = $1;
    }
    if (/PLOIDY\s+\=\s+(\S+)/)
    {
		$controlArgs{PLOIDY} = $1;
    }
    if (/FQ\s+\=\s+(\S+)/)
    {
        $controlArgs{FQ} = $1;
    }
    if (/GENOTYPE_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{GENOTYPE_OUT} = $1;
		if ($passedArgs{runmode} ne "alleles" && $passedArgs{runmode} ne "estPloidy" && $passedArgs{runmode} ne "population2")
		{
			system "mkdir $controlArgs{GENOTYPE_OUT}";
		}
    }
    if (/PHASE_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{PHASE_OUT} = $1;
		if ($passedArgs{runmode} ne "alleles" && $passedArgs{runmode} ne "estPloidy" && $passedArgs{runmode} ne "population2")
		{
			system "mkdir $controlArgs{PHASE_OUT}";
		}
    }
    if (/IUPAC_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{IUPAC_OUT} = $1;
		if ($passedArgs{runmode} ne "alleles" && $passedArgs{runmode} ne "estPloidy" && $passedArgs{runmode} ne "population2")
		{
			system "mkdir $controlArgs{IUPAC_OUT}";
		}
    }
    if (/FASTA_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{FASTA_OUT} = $1;
		if ($passedArgs{runmode} eq "alleles")
		{
			system "mkdir $controlArgs{FASTA_OUT}";
			system "mkdir $controlArgs{FASTA_OUT}/PHASED";
			system "mkdir $controlArgs{FASTA_OUT}/GENOTYPE";
			system "mkdir $controlArgs{FASTA_OUT}/PICKONE";
		}
    }
    if (/SUMMARYSTATS_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{SUMMARYSTATS_OUT} = $1;
		if ($passedArgs{runmode} eq "alleles")
		{
			system "mkdir $controlArgs{SUMMARYSTATS_OUT}";
		}
    }
    if (/ESTPLOIDY_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{ESTPLOIDY_OUT} = $1;
		if ($passedArgs{runmode} eq "estPloidy")
		{
			system "mkdir $controlArgs{ESTPLOIDY_OUT}";
		}
    }
    if (/HELPERSCRIPTS\s+\=\s+(\S+)/)
    {
	    $controlArgs{HELPERSCRIPTS} = $1;
	}
    if (/REF\s+\=\s+(\S+)/)
    {
        $controlArgs{REF} = $1;
    }
    if (/PICARD\s+\=\s+(\S+)/)
    {
        $controlArgs{PICARD} = $1;
        if ($controlArgs{PICARD} =~ m/\S+\.jar/)
        {
        	$controlArgs{PICARD} = "java -jar $controlArgs{PICARD}";
        }	
    }
    if (/BWA\s+\=\s+(\S+)/)
    {
        $controlArgs{BWA} = $1;
    }
    if (/GATK\s+\=\s+(\S+)/)
    {
        $controlArgs{GATK} = $1;
    }
    if (/HPOPG\s+\=\s+(\S+)/)
    {
        $controlArgs{HPOPG} = $1;
    }
    if (/SAMTOOLS\s+\=\s+(\S+)/)
    {
        $controlArgs{SAMTOOLS} = $1;
    }
    if (/BCFTOOLS\s+\=\s+(\S+)/)
    {
        $controlArgs{BCFTOOLS} = $1;
    }
    if (/BAMTOOLS\s+\=\s+(\S+)/)
    {
        $controlArgs{BAMTOOLS} = $1;
    }
    if (/BGZIP\s+\=\s+(\S+)/)
    {
        $controlArgs{BGZIP} = $1;
    }
    if (/TABIX\s+\=\s+(\S+)/)
    {
        $controlArgs{TABIX} = $1;
    }
    if (/SCHEDULER\s+\=\s+(\S+)/)
    {
        $controlArgs{SCHEDULER} = $1;
    }
    if (/GATK_FILTER_EXPRESSION\s+=\s+\"(.+)\"\s+\"(.+)\"/)
    {
    	my $filter = $1;
    	my $filterName = $2;
        push @{$variantFilters{FILTER}}, $filter;
        push @{$variantFilters{FILTERNAME}}, $filterName;
    }
    if (/REFERENCEIND\s+\=\s+(\S+)/)
    {
        $controlArgs{REFERENCEIND} = $1;
    }
    if (/UNIQUE_ONLY\s+\=\s+(\S+)/)
    {
		$controlArgs{UNIQUE_ONLY} = $1;
    }
    if (/REMOVE_READ_DUPLICATES\s+\=\s+(\S+)/)
    {
		$controlArgs{REMOVE_READ_DUPLICATES} = $1;
    }
    if (/OUTPUT_EXPECTED_DOSAGE\s+\=\s+(\S+)/)
    {
		$controlArgs{OUTPUT_EXPECTED_DOSAGE} = $1;
    }
}
close FH1;

if (! exists $variantFilters{FILTER})
{
	print "\n########\nWARNING - no variant filter options specified\n########\n";
	print "The following filters and filter names have been applied by default:\n";
	print "GATK_FILTER_EXPRESSION = \"QD < 2.0\" \"QD_lt2\"\n";
	print "GATK_FILTER_EXPRESSION = \"FS > 60.0\" \"FS_gt60\"\n";
	print "GATK_FILTER_EXPRESSION = \"MQ < 40.0\" \"MQ_lt40\"\n";
	print "GATK_FILTER_EXPRESSION = \"ReadPosRankSum < -8.0\" \"ReadPosRankSum_ltm8\"\n";
	print "GATK_FILTER_EXPRESSION = \"AF < 0.05\" \"AF_lt05\"\n";
	print "GATK_FILTER_EXPRESSION = \"AF > 0.95\" \"AF_gt95\"\n";
	print "GATK_FILTER_EXPRESSION = \"DP < 10\" \"DP_lt10\"\n";
	my @emergencyFilters = ("QD < 2.0","FS > 60.0","MQ < 40.0","ReadPosRankSum < -8.0","AF < 0.05","AF > 0.95","DP < 10");
	my @emergencyFiltersNames = ("QD_lt2","FS_gt60","MQ_lt40","ReadPosRankSum_ltm8","AF_lt05","AF_gt95","DP_lt10");
	for my $i (0..(scalar(@emergencyFilters)-1))
	{
		push @{$variantFilters{FILTER}}, $emergencyFilters[$i];
        push @{$variantFilters{FILTERNAME}}, $emergencyFiltersNames[$i];
	}
}

my %taxaPloidy = ();
open FH1,'<',"$controlArgs{PLOIDY}";
while (<FH1>)
{
    my $line = $_;
    chomp $line;
    my @temp = ();
    @temp = split(/\s+/,$line);
    if (! exists $taxaPloidy{$temp[0]})
    {
		$taxaPloidy{$temp[0]} = $temp[1];
    }
    elsif (exists $taxaPloidy{$temp[0]})
    {   
		die "Duplicate Individual ID $temp[0] in PLOIDY\nThis should be fixed because output data will be lost or incorrect";
    }
}
close FH1;

#altSeqs are used when for --genotype 1
#Otherwise refSeqs are the default when building final fasta files
my %altSeqs = ();
my %refSeqs = ();
my %locusList = ();
my @referenceFasta = glob("$controlArgs{REF}/*.fasta");
foreach my $rf (@referenceFasta)
{
    my $locus = "";
    my $tax = "";
    if ($rf =~ m/$controlArgs{REF}\/(\S+)\.fasta/)
    {
        $locus = $1;
        open FH1,'<',"$rf";
        while (<FH1>)
        {
            if (/^>(\S+)/)
            {
                $tax = $1;
                if ($tax =~ m/(\S+)\-L\d+$/)
                {
                	$tax = $1;
                }
                $refSeqs{$tax}{$locus} = "";
                if (! exists $locusList{$locus})
                {
                    $locusList{$locus} = 1;
                }
            }
            elsif (/(\S+)/)
            {
                my $seq = $1;
                $refSeqs{$tax}{$locus} = $refSeqs{$tax}{$locus} . $seq;
            }
        }
        close FH1;
    }
}

foreach my $tax (sort keys %taxaPloidy)
{
	my @fq1 = ();
	my @fq2 = ();

    if ($passedArgs{runmode} eq "species")
    {
		@fq1 = glob("$controlArgs{FQ}/$tax.R1.*");
		@fq2 = glob("$controlArgs{FQ}/$tax.R2.*");
		system "mkdir $controlArgs{GENOTYPE_OUT}/$tax";
		system "mkdir $controlArgs{PHASE_OUT}/$tax";
		system "mkdir $controlArgs{IUPAC_OUT}/$tax";
		
		open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta";
		foreach my $locus (sort keys %locusList)
		{
	    	if (exists $refSeqs{$tax}{$locus})
		    {
				print OUT1 ">$locus\n$refSeqs{$tax}{$locus}\n";
	    	}
		}
		close OUT1;
		open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.genotype.sh";
		open FH1,'<',"$passedArgs{template}";
		while(<FH1>)
		{
	    	my $line = $_;
	    	chomp $line;
		    $line =~ s/__RUNID__/$tax.id/;
		    $line =~ s/__LOGFILE__/$tax.log/;
	    	print OUT1 "$line\n";
		}
		close FH1;
		print OUT1 "cd $controlArgs{GENOTYPE_OUT}/$tax\n";
	
#####STEP ZERO: Make Reference Databases
#		print "STEP 0: Preparing Reference Databases\n";
		print OUT1 "$controlArgs{PICARD} CreateSequenceDictionary R=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n"; 
		print OUT1 "$controlArgs{BWA} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
		print OUT1 "$controlArgs{SAMTOOLS} faidx $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
	
#####STEP ONE: Map reads
#		print "STEP 1: Mapping Reads\n";
		print OUT1 "$controlArgs{BWA} mem $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta $fq1[0] $fq2[0] | $controlArgs{SAMTOOLS} view -bS -F 4 - | $controlArgs{SAMTOOLS} sort - -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam\n";
	
		print OUT1 "$controlArgs{PICARD} FastqToSam F1=$fq1[0] F2=$fq2[0] O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam SM=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
		print OUT1 "$controlArgs{PICARD} MergeBamAlignment ALIGNED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam UNMAPPED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam R=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
	
#####STEP TWO: Mark duplicates if needed or skip
		if ($controlArgs{REMOVE_READ_DUPLICATES} == 1)
		{
#			print "STEP 2: Marking Duplicates\n";
			print OUT1 "$controlArgs{PICARD} MarkDuplicates I=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam M=$controlArgs{GENOTYPE_OUT}/$tax/$tax.metrics.txt\n";
	
#######STEP THREE: Identify variants, select only SNPs
#			print "STEP 3: Identifying variants\n";
			print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam\n";
			print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -I $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam -ploidy $taxaPloidy{$tax} -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf\n";
		}
		elsif ($controlArgs{REMOVE_READ_DUPLICATES} == 0)
		{
			print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam\n";
			print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -I $controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam -ploidy $taxaPloidy{$tax} -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf\n";
		}
		print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf\n"; 
		print OUT1 "$controlArgs{GATK} VariantAnnotator -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf -A AlleleFraction -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.ab.vcf\n";
		my $gatkfilterstring = "$controlArgs{GATK} VariantFiltration -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.ab.vcf -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.abf.vcf \\";
		for my $fn (0..(scalar(@{$variantFilters{FILTER}})-1))
		{
			if ($fn < (scalar(@{$variantFilters{FILTER}})-1))
			{
				$gatkfilterstring = $gatkfilterstring . "\n --filter-expression \"$variantFilters{FILTER}[$fn]\" --filter-name \"$variantFilters{FILTERNAME}[$fn]\" \\";
			}
			elsif ($fn == (scalar(@{$variantFilters{FILTER}})-1))
			{
				$gatkfilterstring = $gatkfilterstring . "\n --filter-expression \"$variantFilters{FILTER}[$fn]\" --filter-name \"$variantFilters{FILTERNAME}[$fn]\"";
			}
		}
		print OUT1 "$gatkfilterstring\n";
		print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.abf.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf\n"; 

######STEP FOUR: Output new supercontig FASTA with ambiguity codes
#		print "STEP 4: Generating IUPAC FASTA file\n";
		print OUT1 "$controlArgs{GATK} FastaAlternateReferenceMaker -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -O $controlArgs{IUPAC_OUT}/$tax/$tax.iupac.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf --use-iupac-sample $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";

######STEP FIVE: Split BAM and VCF by locus then phase each locus
#		print "STEP 5: Splitting BAM and VCF file by locus\n";
		print OUT1 "$controlArgs{BAMTOOLS} split -in $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam -reference\n";
		print OUT1 "$controlArgs{BGZIP} -c $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf > $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz\n";
		print OUT1 "$controlArgs{TABIX} -p vcf $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz\n";
	
		foreach my $locus (sort keys %locusList)
		{
	    	if (exists $refSeqs{$tax}{$locus})
		    {
				print OUT1 "$controlArgs{TABIX} $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz $locus -h > $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf\n";
				my $locusBAM = "$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.REF_" . "$locus.bam";
				print OUT1 "java -jar $controlArgs{HPOPG} -b $locusBAM -v $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf -p $taxaPloidy{$tax} -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.out -d $controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.log\n";
			}
		}
######STEP SIX: BREAK and distribute on cluster
#		print "STEP 6: Take break to distribute jobs\n";
		system "$controlArgs{SCHEDULER} $controlArgs{GENOTYPE_OUT}/$tax/$tax.genotype.sh";
    }
    
    elsif ($passedArgs{runmode} eq "estPloidy")
    {
		@fq1 = glob("$controlArgs{FQ}/$tax.R1.*");
		@fq2 = glob("$controlArgs{FQ}/$tax.R2.*");
		system "mkdir $controlArgs{ESTPLOIDY_OUT}/$tax";
		
		open OUT1,'>',"$controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta";
		foreach my $locus (sort keys %locusList)
		{
	    	if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		    {
				print OUT1 ">$locus\n$refSeqs{$controlArgs{REFERENCEIND}}{$locus}\n";
	    	}
		}
		close OUT1;
		open OUT1,'>',"$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.estPloidy.sh";
		open FH1,'<',"$passedArgs{template}";
		while(<FH1>)
		{
	    	my $line = $_;
	    	chomp $line;
		    $line =~ s/__RUNID__/$tax.id/;
		    $line =~ s/__LOGFILE__/$tax.log/;
	    	print OUT1 "$line\n";
		}
		close FH1;
		print OUT1 "cd $controlArgs{ESTPLOIDY_OUT}/$tax\n";
	
#####STEP ZERO: Make Reference Databases
#		print "STEP 0: Preparing Reference Databases\n";
		print OUT1 "$controlArgs{PICARD} CreateSequenceDictionary R=$controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n"; 
		print OUT1 "$controlArgs{BWA} index $controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
		print OUT1 "$controlArgs{SAMTOOLS} faidx $controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
	
#####STEP ONE: Map reads
#		print "STEP 1: Mapping Reads\n";
		print OUT1 "$controlArgs{BWA} mem $controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta $fq1[0] $fq2[0] | $controlArgs{SAMTOOLS} view -bS -F 4 - | $controlArgs{SAMTOOLS} sort - -o $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.sorted.bam\n";
	
		print OUT1 "$controlArgs{PICARD} FastqToSam F1=$fq1[0] F2=$fq2[0] O=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.unmapped.bam SM=$controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
		print OUT1 "$controlArgs{PICARD} MergeBamAlignment ALIGNED=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.sorted.bam UNMAPPED=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.unmapped.bam O=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.merged.bam R=$controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
	
#####STEP TWO: Mark duplicates if needed
		if ($controlArgs{REMOVE_READ_DUPLICATES} == 1)
		{
#			print "STEP 2: Marking Duplicates\n";
			print OUT1 "$controlArgs{PICARD} MarkDuplicates I=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.merged.bam O=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.marked.bam M=$controlArgs{ESTPLOIDY_OUT}/$tax/$tax.metrics.txt\n";
	
#######STEP THREE: Identify variants,. Joint genotyping and filtering happens in later step
#			print "STEP 3: Identifying variants\n";
			print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.marked.bam\n";
			print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta -I $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.marked.bam -ploidy 2 -O $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.g.vcf -ERC GVCF\n";
		}
		elsif ($controlArgs{REMOVE_READ_DUPLICATES} == 0)
		{
			print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.merged.bam\n";
			print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{ESTPLOIDY_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta -I $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.merged.bam -ploidy 2 -O $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.g.vcf -ERC GVCF\n";
		}
		
		system "$controlArgs{SCHEDULER} $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.estPloidy.sh";
    }
    
    elsif ($passedArgs{runmode} eq "population1")
    {
		@fq1 = glob("$controlArgs{FQ}/$tax.R1.*");
		@fq2 = glob("$controlArgs{FQ}/$tax.R2.*");
		system "mkdir $controlArgs{GENOTYPE_OUT}/$tax";
		system "mkdir $controlArgs{PHASE_OUT}/$tax";
		system "mkdir $controlArgs{IUPAC_OUT}/$tax";
		
		open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta";
		foreach my $locus (sort keys %locusList)
		{
	    	if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		    {
				print OUT1 ">$locus\n$refSeqs{$controlArgs{REFERENCEIND}}{$locus}\n";
	    	}
		}
		close OUT1;
		open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.pop1.sh";
		open FH1,'<',"$passedArgs{template}";
		while(<FH1>)
		{
	    	my $line = $_;
	    	chomp $line;
		    $line =~ s/__RUNID__/$tax.pop1.id/;
		    $line =~ s/__LOGFILE__/$tax.pop1.log/;
	    	print OUT1 "$line\n";
		}
		close FH1;
		print OUT1 "cd $controlArgs{GENOTYPE_OUT}/$tax\n";
	
#####STEP ZERO: Make Reference Databases
#		print "STEP 0: Preparing Reference Databases\n";
		print OUT1 "$controlArgs{PICARD} CreateSequenceDictionary R=$controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n"; 
		print OUT1 "$controlArgs{BWA} index $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
		print OUT1 "$controlArgs{SAMTOOLS} faidx $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
	
#####STEP ONE: Map reads
#		print "STEP 1: Mapping Reads\n";
		print OUT1 "$controlArgs{BWA} mem $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta $fq1[0] $fq2[0] | $controlArgs{SAMTOOLS} view -bS -F 4 - | $controlArgs{SAMTOOLS} sort - -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam\n";
	
		print OUT1 "$controlArgs{PICARD} FastqToSam F1=$fq1[0] F2=$fq2[0] O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam SM=$controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
		print OUT1 "$controlArgs{PICARD} MergeBamAlignment ALIGNED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam UNMAPPED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam R=$controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta\n";
	
#####STEP TWO: Mark duplicates if needed
		if ($controlArgs{REMOVE_READ_DUPLICATES} == 1)
		{
#			print "STEP 2: Marking Duplicates\n";
			print OUT1 "$controlArgs{PICARD} MarkDuplicates I=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam M=$controlArgs{GENOTYPE_OUT}/$tax/$tax.metrics.txt\n";
	
#######STEP THREE: Identify variants,. Joint genotyping and filtering happens in later step
#			print "STEP 3: Identifying variants\n";
			print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam\n";
			print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta -I $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam -ploidy $taxaPloidy{$tax}  -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.g.vcf -ERC GVCF\n";
		}

		elsif ($controlArgs{REMOVE_READ_DUPLICATES} == 0)
		{	
			print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam\n";
			print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta -I $controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam -ploidy $taxaPloidy{$tax}  -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.g.vcf -ERC GVCF\n";
		}
		system "$controlArgs{SCHEDULER} $controlArgs{GENOTYPE_OUT}/$tax/$tax.pop1.sh";
    }
    
    elsif ($passedArgs{runmode} eq "population2")
    {
    	open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.pop2.sh";
		open FH1,'<',"$passedArgs{template}";
		while(<FH1>)
		{
	    	my $line = $_;
	    	chomp $line;
		    $line =~ s/__RUNID__/$tax.pop2.id/;
		    $line =~ s/__LOGFILE__/$tax.pop2.log/;
	    	print OUT1 "$line\n";
		}
		close FH1;
		print OUT1 "cd $controlArgs{GENOTYPE_OUT}/$tax\n";
    	#######STEP FOUR - continued from pop1: Split the multi-sample vcf into the individual vcf files for processing by the alleles runmode
		print OUT1 "$controlArgs{BCFTOOLS} view -c 1 -O v -s $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf $controlArgs{GENOTYPE_OUT}/joint_genotypes.filtered.biallelic.vcf\n";
		
		######STEP FIVE: Output new supercontig FASTA with ambiguity codes from the split vcf GIVEN the reference sequence
		print OUT1 "$controlArgs{GATK} FastaAlternateReferenceMaker -R $controlArgs{GENOTYPE_OUT}/$tax/$controlArgs{REFERENCEIND}.ref.fasta -O $controlArgs{IUPAC_OUT}/$tax/$tax.iupac.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf\n";
		
		######STEP SIX: Split BAM and VCF by locus then phase each locus
		print OUT1 "$controlArgs{BAMTOOLS} split -in $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam -reference\n";
		print OUT1 "$controlArgs{BGZIP} -c $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf > $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz\n";
		print OUT1 "$controlArgs{TABIX} -p vcf $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz\n";
	
		foreach my $locus (sort keys %locusList)
		{
	    	if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		    {
				print OUT1 "$controlArgs{TABIX} $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz $locus -h > $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf\n";
				my $locusBAM = "$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.REF_" . "$locus.bam";
				print OUT1 "java -jar $controlArgs{HPOPG} -b $locusBAM -v $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf -p $taxaPloidy{$tax} -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.out -d $controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.log\n";
			}
		}
######STEP SIX: BREAK and distribute on cluster
#		print "STEP 6: Take break to distribute jobs\n";
		system "$controlArgs{SCHEDULER} $controlArgs{GENOTYPE_OUT}/$tax/$tax.pop2.sh";
	}
    
    
    elsif ($passedArgs{runmode} eq "alleles")
    {
########
#Do a quick check to see if we are working with population data
		my $ispopulation = 0;
		if (-e "$controlArgs{GENOTYPE_OUT}/intervals.list")
		{
			$ispopulation = 1;
		}
########    	
#    	print "STEP 7: Building phased haplotype sequences for $tax\n";
		if ($passedArgs{genotype} eq "iupac")
		{
		    my $af = "$controlArgs{IUPAC_OUT}/$tax/$tax.iupac.fasta";
		    my $locus = "";
	    	open FH1,'<',"$af";
	    	while (<FH1>)
	    	{
				if (/^>\d+\s+(\S+)\:\d+\-\d+/)
				{
			    	$locus = $1;
			    	$altSeqs{$tax}{$locus} = "";
				}
				elsif (/(\S+)/)
				{
		    		my $seq = $1;
				    $altSeqs{$tax}{$locus} = $altSeqs{$tax}{$locus} . $seq;
				}
	    	}
		}
    	foreach my $locus (sort keys %locusList)
		{
			my $locusdoesexist = 0;
			if ($ispopulation == 0)
			{
		    	if (exists $refSeqs{$tax}{$locus})
		    	{
		    		$locusdoesexist = 1;
		    	}
		    }
		    elsif ($ispopulation == 1)
		    {
		    	if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		    	{
		    		$locusdoesexist = 1;
		    	}
		    }
		    if ($locusdoesexist == 1)
		    {
		    	my $variantsDetected = 0;
		    	open FH1,'<',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf";
		    	while(<FH1>)
		    	{
		    		if (/\S+\s+\d+\s+\.\s+\S+\s+\S+\s+\S+\s+PASS\s+\S+\s+\S+\s+\S+/)
		    		{
		    			$variantsDetected++;
		    		}
		    	}
		    	close FH1;
		    	
		    	if ($variantsDetected > 0)
		    	{
					my $switch = 0;
					my %alleles = ();
					my %posMap = ();
					my $nsnps = 0;
					my %snpList = ();
					my %snpMap = ();
#					print "Processing VCF File - $tax $locus\n";
					open FH1,'<',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf";
					while(<FH1>)
					{
				    	my $line = $_;
					    chomp $line;
			   			if ($switch == 0)
			    		{
							if ($line =~ m/#CHROM\s+.+/)
							{
							    $switch = 1;
							}
		    			}
		    			if ($switch == 1)
		    			{
							if ($line =~ m/\S+\s+(\d+)\s+\.\s+(\S+)\s+(\S+)\s+\S+\s+PASS\s+\S+\s+\S+\s+\S+/)
							{
						    	my $pos = $1;
							    my $a1 = $2;
							    my $a2 = $3;
				    			$alleles{$pos}{0} = $a1;
				    			$alleles{$pos}{1} = $a2;
			    				$nsnps++;
			    				$posMap{$nsnps} = $pos;
			    				$snpList{$pos} = 1;
			    				$snpMap{$pos} = $nsnps;
#			   			 		print "$pos\t$nsnps\t$a1\t$a2\n";
							}
		    			}
					}
					close FH1;
				
					my $fragFile = "$controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.out";
					if (-e $fragFile)
					{
						my %chunks = ();
						my $ploidy = 0;
						my $maxChunk = 0;
						my $nChunks = 0;
						my @snpMatrix = ();
						my %chunkLength = ();
					
						open FH1,'<',"$fragFile";
						while(<FH1>)
						{
			   		 		if (/\#BLOCK\:/)
		    				{
								$nChunks++;
								$chunkLength{$nChunks} = 0;
#								print "Reading Frag - $fragFile BLOCK: $nChunks\n";
					    	}
						    elsif (/^\*.+/)
			    			{
								if ($chunkLength{$nChunks} > $maxChunk)
								{
			    					$maxChunk = $nChunks;
								}
				    		}
					    	elsif (/^(\d+)\s+(.+)/)
			    			{
								my $snp = $1;
								my $vars = $2;
								$chunkLength{$nChunks} = $chunkLength{$nChunks} + 1;
								my @temp = ();
								@temp = split(/\s+/,$vars);
								if ($ploidy == 0)
								{
							    	$ploidy = scalar(@temp);
								}
								push @{$chunks{$nChunks}}, $snp;
								for my $i (0..(scalar(@temp)-1))
								{
							    	$snpMatrix[$snp][$i+1] = $temp[$i];
							    	#print "$snp\t$i\t$temp[$i]\n";
								}
		    				}
						}
						close FH1;


				    	if ($nChunks > 0)
		    			{
	    					open OUT1,'>',"$controlArgs{PHASE_OUT}/$tax/$tax.$locus.phased.fasta";
	    					for my $i (1..$taxaPloidy{$tax})
	    					{
								my $thisSeq = "";
								if ($passedArgs{genotype} eq "consensus")
								{
									if ($ispopulation == 0)
									{
										$thisSeq = $refSeqs{$tax}{$locus};
									}
									elsif ($ispopulation == 1)
									{
										$thisSeq = $refSeqs{$controlArgs{REFERENCEIND}}{$locus};
									}
								}
								elsif ($passedArgs{genotype} eq "iupac")
								{
									$thisSeq = $altSeqs{$tax}{$locus};
								}
								my $phasedSeq = "";
								my @temp = ();
								@temp = split(//,$thisSeq);
								my %varList = ();
								for my $snp (@{$chunks{$maxChunk}})
								{
								    $varList{$snp} = 1;
								}
								for my $j (0..(scalar(@temp)-1))
								{
			    					my $pos = $j+1;
				    				if (! exists $snpMap{$pos})
			    					{
										$phasedSeq = "$phasedSeq" . "$temp[$j]";
		   							}
		    						elsif (exists $snpMap{$pos})
		    						{
										if (exists $varList{$snpMap{$pos}})
										{
#			    							print "Phased SNP Found - $pos c=$temp[$j] - 0=$alleles{$pos}{0} 1=$alleles{$pos}{1}\n";
								    		if ($snpMatrix[$snpMap{$pos}][$i] =~ m/\d+/)
								    		{ 
								    			if ($snpMatrix[$snpMap{$pos}][$i] == 0)                                                          
				    							{	                                                                                                                    
													$phasedSeq = "$phasedSeq" . "$alleles{$pos}{0}";
						    					}                                                                                                                         
											    elsif ($snpMatrix[$snpMap{$pos}][$i] == 1)  
						    					{
													$phasedSeq = "$phasedSeq" . "$alleles{$pos}{1}";
			    								}
											}
											elsif ($snpMatrix[$snpMap{$pos}][$i] !~ m/\d+/)
											{
												if ($snpMatrix[$snpMap{$pos}][$i] eq "-")                                                          
			    								{
			    									$phasedSeq = "$phasedSeq" . "N";
				    							}
				    							else
				    							{
				    								die "Critical Error! SNP Matrix from HPOPG has characters --- $snpMatrix[$snpMap{$pos}][$i] --- it should not\n Stopping until problems are resolved\n";
			    								}									
											}
										}
										elsif (! exists $varList{$snpMap{$pos}})
										{
#			    							print "Phased SNP Found on Short end - $pos c=$temp[$j] - 0=$alleles{$pos}{0} 1=$alleles{$pos}{1}\n";
								    		if ($passedArgs{genotype} eq "consensus")
											{
								    			$phasedSeq = "$phasedSeq" . "N";
								    		}
								    		elsif ($passedArgs{genotype} eq "iupac")
											{
								    			$phasedSeq = "$phasedSeq" . "$temp[$j]";
								    		}
										}
			    					}
								}
								my $header = "$tax" . "__" . "$locus" . "__" . "$i";
								print OUT1 ">$header\n$phasedSeq\n";
	    					}
		   					close OUT1;
		    			}
		    		}
		    	}
			}
    	}    	
	}
}


if ($passedArgs{runmode} eq "estPloidy")
{
	open OUT1,'>',"$controlArgs{ESTPLOIDY_OUT}/intervals.list";
	foreach my $locus (sort keys %locusList)
	{
		if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		{
	    	print OUT1 "$locus\n";
	    }
	}
	open OUT1,'>',"$controlArgs{ESTPLOIDY_OUT}/jointGenotyping.estPloidy.sh";
	open FH1,'<',"$passedArgs{template}";
	while(<FH1>)
	{
	   	my $line = $_;
	   	chomp $line;
	    $line =~ s/__RUNID__/jointGenotyping.estPloidy.id/;
	    $line =~ s/__LOGFILE__/jointGenotyping.estPloidy.log/;
	   	print OUT1 "$line\n";
	}
	close FH1;
	print OUT1 "cd $controlArgs{ESTPLOIDY_OUT}\n";
	print OUT1 "$controlArgs{GATK} GenomicsDBImport";
	foreach my $tax (sort keys %taxaPloidy)
	{
		print OUT1 " -V $controlArgs{ESTPLOIDY_OUT}/$tax/$tax.g.vcf";
	}
	print OUT1 " --genomicsdb-workspace-path joint_estPloidy --intervals intervals.list\n";
	print OUT1 "$controlArgs{GATK} GenotypeGVCFs -R $controlArgs{ESTPLOIDY_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V gendb://joint_estPloidy -O $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.vcf\n";
	print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{ESTPLOIDY_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.biallelic.vcf\n"; 
	my $gatkfilterstring = "$controlArgs{GATK} VariantFiltration -R $controlArgs{ESTPLOIDY_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.biallelic.vcf -O $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.filtered.vcf \\";
	for my $fn (0..(scalar(@{$variantFilters{FILTER}})-1))
	{
		if ($fn < (scalar(@{$variantFilters{FILTER}})-1))
		{
			$gatkfilterstring = $gatkfilterstring . "\n --filter-expression \"$variantFilters{FILTER}[$fn]\" --filter-name \"$variantFilters{FILTERNAME}[$fn]\" \\";
		}
		elsif ($fn == (scalar(@{$variantFilters{FILTER}})-1))
		{
			$gatkfilterstring = $gatkfilterstring . "\n --filter-expression \"$variantFilters{FILTER}[$fn]\" --filter-name \"$variantFilters{FILTERNAME}[$fn]\"";
		}
	}
	print OUT1 "$gatkfilterstring\n";
	print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{ESTPLOIDY_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.filtered.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.filtered.biallelic.vcf\n"; 

#######STEP FOUR: Collect ratio of alt reads to ref reads and estimate ploidy with normal mixture model
	my $tempsched = "scheduler_" . "$controlArgs{SCHEDULER}";
	print OUT1 "perl $controlArgs{HELPERSCRIPTS}/getAB.pl $controlArgs{ESTPLOIDY_OUT}/joint_estPloidy.filtered.biallelic.vcf $controlArgs{ESTPLOIDY_OUT}\n";
	foreach my $tax (sort keys %taxaPloidy)
    {
	    print OUT1 "perl $controlArgs{HELPERSCRIPTS}/estimatePloidy.pl $tempsched $controlArgs{PHASING_ROOT}/$passedArgs{template} $controlArgs{ESTPLOIDY_OUT} $controlArgs{HELPERSCRIPTS} $tax\n";
	}
}

if ($passedArgs{runmode} eq "population1")
{
	open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/intervals.list";
	foreach my $locus (sort keys %locusList)
	{
		if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		{
	    	print OUT1 "$locus\n";
	    }
	}
	open OUT1,'>',"$controlArgs{GENOTYPE_OUT}/jointGenotyping.sh";
	open FH1,'<',"$passedArgs{template}";
	while(<FH1>)
	{
	   	my $line = $_;
	   	chomp $line;
	    $line =~ s/__RUNID__/jointGenotyping.id/;
	    $line =~ s/__LOGFILE__/jointGenotyping.log/;
	   	print OUT1 "$line\n";
	}
	close FH1;
	print OUT1 "cd $controlArgs{GENOTYPE_OUT}\n";
	print OUT1 "$controlArgs{GATK} GenomicsDBImport";
	foreach my $tax (sort keys %taxaPloidy)
	{
		print OUT1 " -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.g.vcf";
	}
	print OUT1 " --genomicsdb-workspace-path joint_genotypes --intervals intervals.list\n";
	print OUT1 "$controlArgs{GATK} GenotypeGVCFs -R $controlArgs{GENOTYPE_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V gendb://joint_genotypes -O $controlArgs{GENOTYPE_OUT}/joint_genotypes.vcf\n";
	print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V $controlArgs{GENOTYPE_OUT}/joint_genotypes.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/joint_genotypes.biallelic.vcf\n"; 
	my $gatkfilterstring = "$controlArgs{GATK} VariantFiltration -R $controlArgs{GENOTYPE_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V $controlArgs{GENOTYPE_OUT}/joint_genotypes.biallelic.vcf -O $controlArgs{GENOTYPE_OUT}/joint_genotypes.filtered.vcf \\";
	for my $fn (0..(scalar(@{$variantFilters{FILTER}})-1))
	{
		if ($fn < (scalar(@{$variantFilters{FILTER}})-1))
		{
			$gatkfilterstring = $gatkfilterstring . "\n --filter-expression \"$variantFilters{FILTER}[$fn]\" --filter-name \"$variantFilters{FILTERNAME}[$fn]\" \\";
		}
		elsif ($fn == (scalar(@{$variantFilters{FILTER}})-1))
		{
			$gatkfilterstring = $gatkfilterstring . "\n --filter-expression \"$variantFilters{FILTER}[$fn]\" --filter-name \"$variantFilters{FILTERNAME}[$fn]\"";
		}
	}
	print OUT1 "$gatkfilterstring\n";
	print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$controlArgs{REFERENCEIND}/$controlArgs{REFERENCEIND}.ref.fasta -V $controlArgs{GENOTYPE_OUT}/joint_genotypes.filtered.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/joint_genotypes.filtered.biallelic.vcf\n"; 

	#This job is not automatically submitted since it needs to wait for the individual genotype outputs
	print "\nWhen all individual genotyping and gvcf files are completed, submit the joint genotyping file $controlArgs{GENOTYPE_OUT}/jointGenotyping.sh\nYou may consider adjusting the RAM and cpu configuration for this one.\n";
}

if ($passedArgs{runmode} eq "alleles")
{
#    print "STEP 8: Output fasta files of orthologous phased haplotype sequences\n";
	my %phasedSeq = ();
	my %tax2phase = ();
	
	foreach my $tax (sort keys %taxaPloidy)
	{
    	my @phasedFastaFiles = ();
    	@phasedFastaFiles = glob("$controlArgs{PHASE_OUT}/$tax/$tax.*.phased.fasta");
    	foreach my $pff (@phasedFastaFiles)
    	{
			if ($pff =~ m/$controlArgs{PHASE_OUT}\/$tax\/$tax\.(\S+)\.phased\.fasta/)
			{
	    		my $locus = $1;
			    my %tempSeq = ();
			    my $allele = "";
	    		open FH1,'<',"$pff";
	    		while (<FH1>)
	   			{
					if (/^>(\S+)/)
					{
					    $allele = $1;
		    			$allele =~ s/__\S+__/__/;
					}
					elsif (/(\S+)/)
					{
		    			my $seq = $1;
		   	 			my $switch = 0;
					    foreach $allele (keys %tempSeq)
		    			{
							if ($seq eq $tempSeq{$allele})
							{
########
#The UNIQUE_ONLY flag works here by keeping the switch variable always off
								if ($controlArgs{UNIQUE_ONLY} == 1)
								{
									$switch = 1;
								}
########
							}
		    			}
		    			if ($switch == 0)
		    			{
							$tempSeq{$allele} = $seq;
							$phasedSeq{$locus}{$allele} = $seq;
							push @{$tax2phase{$locus}{$tax}}, $allele;
		   				}
					}
	    		}
			    close FH1;
			}
    	}
	}

########
#Added the pick-one fasta output files 20220214
#Decided to let randomly chosen alleles retain numbered headers and somebody can always choose to remove the __suffix
########

	foreach my $locus (sort keys %locusList)
	{
    	open OUT1,'>',"$controlArgs{FASTA_OUT}/PHASED/$locus.fasta";
    	open OUT2,'>',"$controlArgs{FASTA_OUT}/PICKONE/$locus.fasta";
    	foreach my $tax (sort keys %taxaPloidy)
    	{
    		if ($passedArgs{genotype} eq "consensus")
    		{
				if (exists $refSeqs{$tax}{$locus})
				{
				    if (exists $tax2phase{$locus}{$tax}[0])
	    			{
						for my $i (0..(scalar(@{$tax2phase{$locus}{$tax}}) - 1))
						{
					    	print OUT1 ">$tax2phase{$locus}{$tax}[$i]\n$phasedSeq{$locus}{$tax2phase{$locus}{$tax}[$i]}\n";
						}
						my $randomAllele = int(rand(scalar(@{$tax2phase{$locus}{$tax}})));
						print OUT2 ">$tax2phase{$locus}{$tax}[$randomAllele]\n$phasedSeq{$locus}{$tax2phase{$locus}{$tax}[$randomAllele]}\n";
		    		}
	    			else
				    {
				    	if ($controlArgs{OUTPUT_EXPECTED_DOSAGE} == 0)
				    	{
				    		my $referenceHeader = "$tax" . "__REF";
							print OUT1 ">$referenceHeader\n$refSeqs{$tax}{$locus}\n";
							print OUT2 ">$referenceHeader\n$refSeqs{$tax}{$locus}\n";
						}
						elsif ($controlArgs{OUTPUT_EXPECTED_DOSAGE} == 1)
				    	{
				    		for my $i (1..($taxaPloidy{$tax}))
				    		{
			    				my $header = "$tax" . "__" . "$i";
								print OUT1 ">$header\n$refSeqs{$tax}{$locus}\n";
	    					}
	    					my $referenceHeader = "$tax" . "__REF";
	    					print OUT2 ">$referenceHeader\n$refSeqs{$tax}{$locus}\n";
				    	}
	    			}
				}
			}
			if ($passedArgs{genotype} eq "iupac")
    		{
				if (exists $altSeqs{$tax}{$locus})
				{
				    if (exists $tax2phase{$locus}{$tax}[0])
	    			{
						for my $i (0..(scalar(@{$tax2phase{$locus}{$tax}}) - 1))
						{
					    	print OUT1 ">$tax2phase{$locus}{$tax}[$i]\n$phasedSeq{$locus}{$tax2phase{$locus}{$tax}[$i]}\n";
						}
						my $randomAllele = int(rand(scalar(@{$tax2phase{$locus}{$tax}})));
						print OUT2 ">$tax2phase{$locus}{$tax}[$randomAllele]\n$phasedSeq{$locus}{$tax2phase{$locus}{$tax}[$randomAllele]}\n";

		    		}
	    			else
				    {
				    	if ($controlArgs{OUTPUT_EXPECTED_DOSAGE} == 0)
				    	{
				    		my $referenceHeader = "$tax" . "__REF";
							print OUT1 ">$referenceHeader\n$altSeqs{$tax}{$locus}\n";
							print OUT2 ">$referenceHeader\n$altSeqs{$tax}{$locus}\n";
	    				}
	    				elsif ($controlArgs{OUTPUT_EXPECTED_DOSAGE} == 1)
	    				{
	    					for my $i (1..($taxaPloidy{$tax}))
				    		{
			    				my $header = "$tax" . "__" . "$i";
								print OUT1 ">$header\n$altSeqs{$tax}{$locus}\n";
								print OUT2 ">$header\n$altSeqs{$tax}{$locus}\n";
	    					}
	    				}
	    			}
				}
			}
    	}
    	close OUT1;
    	close OUT2;
	}
	
########
#Added per-locus genotype fasta files to FASTA_OUT 20201109
########	
#	print "STEP 9: Output fasta files of orthologous unphased genotype sequences\n";
	my %genotypeSeq = ();
	foreach my $tax (sort keys %taxaPloidy)
	{
		open FH1,'<',"$controlArgs{IUPAC_OUT}/$tax/$tax.iupac.fasta";
		my $locus = "";
		while(<FH1>)
		{
			if (/^>\S+\s+(\S+)\:\d+\-\d+/)
			{
				$locus = $1;
				$genotypeSeq{$tax}{$locus} = ""; 
			}
			elsif (/(\S+)/)
			{
				my $seq = $1;
				$genotypeSeq{$tax}{$locus} = $genotypeSeq{$tax}{$locus} . $seq;
			}
		}
		close FH1;
	}
	
	foreach my $locus (sort keys %locusList)
	{
		open OUT1,'>',"$controlArgs{FASTA_OUT}/GENOTYPE/$locus.fasta";
    	foreach my $tax (sort keys %taxaPloidy)
    	{
    		if (exists $genotypeSeq{$tax}{$locus})
			{
    			print OUT1 ">$tax\n$genotypeSeq{$tax}{$locus}\n";
    		}
    	}
    	close OUT1;
    }
    
########
#Added basic phasing stats to SUMMARYSTATS_OUT
########
#	print "STEP 10: Collecting Phasing Statistics\n";
	my %phasein = ();
	my %alleleCounts = ();
	my $maxAlleles = 0;
	my %globalStats = ();
	my %globalCounts = ();

	foreach my $tax (sort keys %taxaPloidy)
	{
		if ($taxaPloidy{$tax} > $maxAlleles)
		{
			$maxAlleles = $taxaPloidy{$tax};	
		}
    	foreach my $locus (sort keys %locusList)
	    {
			my $nvar = 0;
			my $nblocks = 0;
			my $het = 0;
			my $bl = 0;
			my $lbl = 0;
			
################################
#Do a quick check to see if we are working with population data
			my $ispopulation = 0;
			if (-e "$controlArgs{GENOTYPE_OUT}/intervals.list")
			{
				$ispopulation = 1;
			}
################################
################################
#And a slightly more extraneous way to check if the locus exists based on if we are looking at the de novo or reference sequences
			my $locusdoesexist = 0;
			if ($ispopulation == 0)
			{
		    	if (exists $refSeqs{$tax}{$locus})
		    	{
		    		$locusdoesexist = 1;
		    	}
		    }
		    elsif ($ispopulation == 1)
		    {
		    	if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		    	{
		    		$locusdoesexist = 1;
		    	}
		    }
################################
			if ($locusdoesexist == 1)
			{
	    		open FH1,'<',"$controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf";
			    while (<FH1>)
	    		{
					if (/PASS/)
					{
					    $nvar++;
					}
	    		}
			    close FH1;

				my $fragFile = "$controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.out";
				if (-e $fragFile)
				{
					open FH1,'<',"$fragFile";
					#BLOCK: offset: 1 len: 19 phased: 76
					while (<FH1>)
					{
					    if(/#BLOCK:\s+offset:\s+\d+\s+len:\s+(\d+)\s+phased:\s+\d+/)
		    			{
							$bl = $1;
							if ($bl > $lbl)
							{
			    				$lbl = $bl;
							}
							$nblocks++;
		    			}
					}
					close FH1;
					
					if ($ispopulation == 0)
					{
		    			$het = $nvar/(length($refSeqs{$tax}{$locus}));
		    		}
		    
		    		elsif ($ispopulation == 1)
		   	 		{
		    			$het = $nvar/(length($refSeqs{$controlArgs{REFERENCEIND}}{$locus}));
		    		}
		    	}
			}
			$phasein{$tax}{$locus}{nvar} = $nvar;
			$phasein{$tax}{$locus}{het} = $het;
			$phasein{$tax}{$locus}{nblocks} = $nblocks;
			$phasein{$tax}{$locus}{lbl} = $lbl;
		}
	}

	foreach my $locus (sort keys %locusList)
	{
    	open FH1,'<',"$controlArgs{FASTA_OUT}/PHASED/$locus.fasta";
    	while (<FH1>)
    	{
			if (/^>(\S+)\_\_\S+/)
			{
			    my $tax = $1;
	   			if (! exists $alleleCounts{$tax}{$locus})
			    {
					$alleleCounts{$tax}{$locus} = 1;
	   			}
	   			elsif (exists $alleleCounts{$tax}{$locus})
	   			{   
					$alleleCounts{$tax}{$locus} = $alleleCounts{$tax}{$locus} + 1;
	   			}
			}
    	}
    	close FH1;
	}

########
#Print out all per-locus stats per individual as well as the global counts that are average over present loci
########
	
	foreach my $tax (sort keys %taxaPloidy)
	{
    	open OUT1,'>',"$controlArgs{SUMMARYSTATS_OUT}/$tax.phasingSummary.txt";
		print OUT1 "LOCUS\tLENGTH\tNVAR\tHET\tNBLOCKS\tLONGESTBL";
		for my $i (1..$maxAlleles)
		{
			print OUT1 "\t$i";
			$globalCounts{$tax}{$i} = 0;
		}
		print OUT1 "\n";
		$globalStats{$tax}{nloci} = 0;
		$globalStats{$tax}{ninvarloci} = 0;
		$globalStats{$tax}{nvarloci} = 0;
		$globalStats{$tax}{nphaseloci} = 0;
		$globalStats{$tax}{loclen} = 0;
		$globalStats{$tax}{nvar} = 0;
		$globalStats{$tax}{het} = 0;
		$globalStats{$tax}{nblocks} = 0;
		$globalStats{$tax}{lbl} = 0;
		
	    foreach my $locus (sort keys %locusList)
    	{
################################
#Do a quick check to see if we are working with population data
			my $ispopulation = 0;
			if (-e "$controlArgs{GENOTYPE_OUT}/intervals.list")
			{
				$ispopulation = 1;
			}
################################
################################
#And a slightly more extraneous way to check if the locus exists based on if we are looking at the de novo or reference sequences
			my $locusdoesexist = 0;
			if ($ispopulation == 0)
			{
		    	if (exists $refSeqs{$tax}{$locus})
		    	{
		    		$locusdoesexist = 1;
		    	}
		    }
		    elsif ($ispopulation == 1)
		    {
		    	if (exists $refSeqs{$controlArgs{REFERENCEIND}}{$locus})
		    	{
		    		($locusdoesexist = 1);
		    	}
		    }
################################
			if ($locusdoesexist == 1)
			{	
				my $refLen = 0;
				if ($ispopulation == 0)
				{
					$refLen = length($refSeqs{$tax}{$locus});
				}
				elsif ($ispopulation == 1)
				{
					$refLen = length($refSeqs{$controlArgs{REFERENCEIND}}{$locus});
				}
	    		
	    		if (exists $alleleCounts{$tax}{$locus})
	    		{
	    			print OUT1 "$locus\t$refLen\t$phasein{$tax}{$locus}{nvar}\t$phasein{$tax}{$locus}{het}\t$phasein{$tax}{$locus}{nblocks}\t$phasein{$tax}{$locus}{lbl}";
	    			for my $i (1..$maxAlleles)
	    			{
	    				if ($alleleCounts{$tax}{$locus} == $i)
	    				{
							print OUT1 "\t1";
							$globalCounts{$tax}{$i} = $globalCounts{$tax}{$i} + 1;
						}
						else
						{
							print OUT1 "\t0";
						}
	    			}
					print OUT1 "\n";
				
				
					$globalStats{$tax}{nloci} = $globalStats{$tax}{nloci} + 1;
					$globalStats{$tax}{loclen} = $globalStats{$tax}{loclen} + $refLen;
				
					if ($phasein{$tax}{$locus}{nblocks} > 0 && $phasein{$tax}{$locus}{nvar} > 0)
					{
						$globalStats{$tax}{nvar} = $globalStats{$tax}{nvar} + $phasein{$tax}{$locus}{nvar};
						$globalStats{$tax}{het} = $globalStats{$tax}{het} + $phasein{$tax}{$locus}{het};
						$globalStats{$tax}{nblocks} = $globalStats{$tax}{nblocks} + $phasein{$tax}{$locus}{nblocks};
						$globalStats{$tax}{lbl} = $globalStats{$tax}{lbl} + $phasein{$tax}{$locus}{lbl};
						$globalStats{$tax}{nphaseloci} = $globalStats{$tax}{nphaseloci} + 1;
						$globalStats{$tax}{nvarloci} = $globalStats{$tax}{nvarloci} + 1;
					}
					if ($phasein{$tax}{$locus}{nblocks} == 0 && $phasein{$tax}{$locus}{nvar} == 0)
					{
						$globalStats{$tax}{ninvarloci} = $globalStats{$tax}{ninvarloci} + 1;
					}
					if ($phasein{$tax}{$locus}{nblocks} == 0 && $phasein{$tax}{$locus}{nvar} > 0)
					{
						$globalStats{$tax}{nvar} = $globalStats{$tax}{nvar} + $phasein{$tax}{$locus}{nvar};
						$globalStats{$tax}{het} = $globalStats{$tax}{het} + $phasein{$tax}{$locus}{het};
						$globalStats{$tax}{nvarloci} = $globalStats{$tax}{nvarloci} + 1;
					}
					if ($phasein{$tax}{$locus}{nblocks} > 0 && $phasein{$tax}{$locus}{nvar} == 0)
					{
						print "Critical error: Cannot have phased loci without passing variants!\n"
					}
				}
				if (! exists $alleleCounts{$tax}{$locus})
	    		{
	    			print OUT1 "$locus\t0\t0\t0\t0\t0";
					for my $i (1..$maxAlleles)
	    			{
	    				print OUT1 "\tNA";
		    		}
		    		print OUT1 "\n";
	    		}
			}
			elsif ($locusdoesexist == 0)
			{
				print OUT1 "$locus\t0\t0\t0\t0\t0";
				for my $i (1..$maxAlleles)
	    		{
	    			print OUT1 "\tNA";
	    		}
	    		print OUT1 "\n";
			}
    	}
    	close OUT1;
    }
    
    open OUT1,'>',"$controlArgs{SUMMARYSTATS_OUT}/averagePhasingStats.txt";
    print OUT1 "INDIVIDUAL\tNLOCI\tNVARLOCI\tNINVLOCI\tNPHASELOCI\tAVG_LENGTH\tAVG_NVAR\tAVG_HET\tAVG_NBLOCKS\tAVG_LONGESTBL";
    for my $i (1..$maxAlleles)
	{
		print OUT1 "\t$i";
	}
	print OUT1 "\n";
	foreach my $tax (sort keys %taxaPloidy)
	{
		if ($globalStats{$tax}{nloci} > 0)
		{
			$globalStats{$tax}{loclen} = $globalStats{$tax}{loclen}/$globalStats{$tax}{nloci};
			$globalStats{$tax}{nvar} = $globalStats{$tax}{nvar}/$globalStats{$tax}{nloci};
			$globalStats{$tax}{het} = $globalStats{$tax}{het}/$globalStats{$tax}{nloci};
			if ($globalStats{$tax}{nphaseloci} > 0)
			{
				$globalStats{$tax}{nblocks} = $globalStats{$tax}{nblocks}/$globalStats{$tax}{nphaseloci};
				$globalStats{$tax}{lbl} = $globalStats{$tax}{lbl}/$globalStats{$tax}{nphaseloci};
			}
			print OUT1 "$tax\t$globalStats{$tax}{nloci}\t$globalStats{$tax}{nvarloci}\t$globalStats{$tax}{ninvarloci}\t$globalStats{$tax}{nphaseloci}\t$globalStats{$tax}{loclen}\t$globalStats{$tax}{nvar}\t$globalStats{$tax}{het}\t$globalStats{$tax}{nblocks}\t$globalStats{$tax}{lbl}";
		}
		elsif ($globalStats{$tax}{nloci} == 0)
		{
			print OUT1 "$tax\t0\t0\t0\t0\t0\t0\t0\t0\t0";
		}
		for my $i (1..$maxAlleles)
		{
			if ($i <= $taxaPloidy{$tax})
			{
				print OUT1 "\t$globalCounts{$tax}{$i}";
			}
			elsif ($i > $taxaPloidy{$tax})
			{
				print OUT1 "\tNA";
			}
		}
		print OUT1 "\n";
	}
	close OUT1;
}
exit;
