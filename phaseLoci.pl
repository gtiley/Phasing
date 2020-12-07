#!/usr/bin/perl -w
use strict;

#----------------------------------------------------------------------------------------#
#George P. Tiley and Andrew A. Crowl
#09 Novermber 2020
#contact: george.tiley@duke.edu
#contact: andrew.crowl@duke.edu
#Genotype target-enrichment loci of known ploidy
#Phase haplotype sequences for each locus
#----------------------------------------------------------------------------------------#

# Accepted AND necessary commands
my @checkArgs = ("controlFile","runmode","template","genotype");
my %passedArgs = ();
if (scalar(@ARGV) == 0)
{
die "/*--------INPUT PARAMETERS--------*/\n
--controlFile STRING <control file>
--runmode INT <0 or 1 or 2>
--template STRING <template file for distributing on a cluster>
--genotype INT <0 or 1>

\n/*--------EXAMPLE COMMAND--------*/\n
    perl phaseLoci.pl --controlFile phase.ctl --runmode 1 --template template.sh --genotype 0\n

\n/*--------FLAG OPTIONS--------*/\n
--runmode
	0 = run each step serially (a string is still needed for --template although it is not used)
	1 = distribute each sample seperatley on a cluster with directives given in template.sh
	2 = Use phasing information to make orthologous fasta files with haplotype sequences AFTER running --bacth as 0 or 1
	
--genotype
	0 = Use the original reference sequence for each sample for each locus when there are no variants to phase
	1 = Use genotyped sequence with IUPAC codes when there are variants, but cannot be phased\n
	
\n/*--------NOTES--------*/\n
This scripts automates a lot of the genotyping and phasing process, but has to be ran twice.
The first time use --runmode 0 OR --runmode 1 depending on if you are running everything serially on a local machine or distributing on a cluster.
All steps assume a single processor is used or allocated. This takes about 30 minutes for a sample of >100x but <1000x coverage. If you want to change these options you may have to alter the appropriate lines relevant to BWA and GATK.
We have yet to explore effects of different hard filtering strategies through GATK, this might require fine-tuning for low (e.g. <=15x) or very high (>1000x) depth data.";

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
		}
		if ($ARGV[$i] eq "--template")
		{
        	$passedArgs{template} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "--genotype")
        {
		    $passedArgs{genotype} = $ARGV[$i+1];
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
open FH1,'<',"$passedArgs{controlFile}";
while (<FH1>)
{
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
		if ($passedArgs{runmode} < 2)
		{
			system "mkdir $controlArgs{GENOTYPE_OUT}";
		}
    }
    if (/PHASE_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{PHASE_OUT} = $1;
		if ($passedArgs{runmode} < 2)
		{
			system "mkdir $controlArgs{PHASE_OUT}";
		}
    }
    if (/IUPAC_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{IUPAC_OUT} = $1;
		if ($passedArgs{runmode} < 2)
		{
			system "mkdir $controlArgs{IUPAC_OUT}";
		}
    }
    if (/FASTA_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{FASTA_OUT} = $1;
		if ($passedArgs{runmode} == 2)
		{
			system "mkdir $controlArgs{FASTA_OUT}";
			system "mkdir $controlArgs{FASTA_OUT}/PHASED";
			system "mkdir $controlArgs{FASTA_OUT}/GENOTYPE";
		}
    }
    if (/SUMMARYSTATS_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{SUMMARYSTATS_OUT} = $1;
		if ($passedArgs{runmode} == 2)
		{
			system "mkdir $controlArgs{SUMMARYSTATS_OUT}";
		}
    }
    if (/REF\s+\=\s+(\S+)/)
    {
        $controlArgs{REF} = $1;
    }
    if (/PICARD\s+\=\s+(\S+)/)
    {
        $controlArgs{PICARD} = $1;
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
}
close FH1;

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
    }
}

foreach my $tax (sort keys %taxaPloidy)
{
	my @fq1 = ();
	my @fq2 = ();
	if ($passedArgs{runmode} == 0)
    {
    	print "####--------Running In Serial Mode--------####\n";
    	print "Warning - Each command will happen subsequently for genotyping - this may take a long time if there are many samples\n";
    	print "Warning - You may be running this locally - Allow ~500MB of disk space for each individual until phasing and clean-up are done\n";
    	
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
	
		print "Genotyping and Phasing Individual - $tax\n";
#####STEP ZERO: Make Reference Databases
		print "STEP 0: Preparing Reference Databases\n";
		system "$controlArgs{PICARD} CreateSequenceDictionary R=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta"; 
		system "$controlArgs{BWA} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta";
		system "$controlArgs{SAMTOOLS} faidx $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta";
	
#####STEP ONE: Map reads
		print "STEP 1: Mapping Reads\n";
		system"$controlArgs{BWA} mem $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta $fq1[0] $fq2[0] | $controlArgs{SAMTOOLS} view -bS -F 4 - | $controlArgs{SAMTOOLS} sort - -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam";
	
		system "$controlArgs{PICARD} FastqToSam F1=$fq1[0] F2=$fq2[0] O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam SM=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta";
		system "$controlArgs{PICARD} MergeBamAlignment ALIGNED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam UNMAPPED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam R=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta";
	
#####STEP TWO: Mark duplicates
		print "STEP 2: Marking Duplicates\n";
		system "$controlArgs{PICARD} MarkDuplicates I=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam M=$controlArgs{GENOTYPE_OUT}/$tax/$tax.metrics.txt";
	
#######STEP THREE: Identify variants, select only SNPs
		print "STEP3: Identifying variants\n";
		system "$controlArgs{SAMTOOLS} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam";
		system "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -I $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam -ploidy $taxaPloidy{$tax} -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf";
		system "$controlArgs{GATK} VariantAnnotator -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf -A AlleleFraction -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.ab.vcf";
		system "$controlArgs{GATK} VariantFiltration -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.ab.vcf -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.abf.vcf \\
                --filter-expression \"QD < 2.0\" --filter-name \"QD_lt2\" \\
                --filter-expression \"FS > 60.0\" --filter-name \"FS_gt60\" \\
                --filter-expression \"MQ < 40.0\" --filter-name \"MQ_lt40\" \\
                --filter-expression \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum_ltm8\" \\
                --filter-expression \"AF < 0.05\" --filter-name \"AF_05\" \\
                --filter-expression \"AF > 0.95\" --filter-name \"AF_05\"";
		system "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.abf.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf"; 

######STEP FOUR: Output new supercontig FASTA with ambiguity codes
		print "STEP 4: Generating IUPAC FASTA file\n";
		system "$controlArgs{GATK} FastaAlternateReferenceMaker -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -O $controlArgs{IUPAC_OUT}/$tax/$tax.iupac.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf --use-iupac-sample $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";

######STEP FIVE: Split BAM and VCF by locus then phase each locus
		print "STEP 5: Splitting BAM and VCF file by locus\n";
		system "$controlArgs{BAMTOOLS} split -in $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam -reference";
		system "$controlArgs{BGZIP} -c $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf > $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz";
		system "$controlArgs{TABIX} -p vcf $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz";
	
		foreach my $locus (sort keys %locusList)
		{
	    	if (exists $refSeqs{$tax}{$locus})
		    {
				system "$controlArgs{TABIX} $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz $locus -h > $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf";
				my $locusBAM = "$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.REF_" . "$locus.bam";
				system "java -jar $controlArgs{HPOPG} -b $locusBAM -v $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.$locus.vcf -p $taxaPloidy{$tax} -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.out -d $controlArgs{GENOTYPE_OUT}/$tax/$tax.$locus.phase.log";
			}
		}
######STEP SIX: Cleanup
		print "STEP 6: Cleaning up unecessary reference files and alignments used for genotyping\n";
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam.bai");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.dict");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta.amb");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta.ann");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta.bwt");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta.fai");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta.pac");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta.sa");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf.idx");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.idx");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.gz");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf.tbi");
		unlink("$controlArgs{GENOTYPE_OUT}/$tax/xieminzhuBAM2frag.txt");
    }

    elsif ($passedArgs{runmode} == 1)
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
		print "STEP 0: Preparing Reference Databases\n";
		print OUT1 "$controlArgs{PICARD} CreateSequenceDictionary R=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n"; 
		print OUT1 "$controlArgs{BWA} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
		print OUT1 "$controlArgs{SAMTOOLS} faidx $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
	
#####STEP ONE: Map reads
		print "STEP 1: Mapping Reads\n";
		print OUT1 "$controlArgs{BWA} mem $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta $fq1[0] $fq2[0] | $controlArgs{SAMTOOLS} view -bS -F 4 - | $controlArgs{SAMTOOLS} sort - -o $controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam\n";
	
		print OUT1 "$controlArgs{PICARD} FastqToSam F1=$fq1[0] F2=$fq2[0] O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam SM=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
		print OUT1 "$controlArgs{PICARD} MergeBamAlignment ALIGNED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.sorted.bam UNMAPPED=$controlArgs{GENOTYPE_OUT}/$tax/$tax.unmapped.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam R=$controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";
	
#####STEP TWO: Mark duplicates
		print "STEP 2: Marking Duplicates\n";
		print OUT1 "$controlArgs{PICARD} MarkDuplicates I=$controlArgs{GENOTYPE_OUT}/$tax/$tax.merged.bam O=$controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam M=$controlArgs{GENOTYPE_OUT}/$tax/$tax.metrics.txt\n";
	
#######STEP THREE: Identify variants, select only SNPs
		print "STEP 3: Identifying variants\n";
		print OUT1 "$controlArgs{SAMTOOLS} index $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam\n";
		print OUT1 "$controlArgs{GATK} HaplotypeCaller -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -I $controlArgs{GENOTYPE_OUT}/$tax/$tax.marked.bam -ploidy $taxaPloidy{$tax} -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf\n";
		print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf\n"; 
		print OUT1 "$controlArgs{GATK} VariantAnnotator -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.vcf -A AlleleFraction -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.ab.vcf\n";
		print OUT1 "$controlArgs{GATK} VariantFiltration -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.ab.vcf -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.abf.vcf \\
                       --filter-expression \"QD < 2.0\" --filter-name \"QD_lt2\" \\
                       --filter-expression \"FS > 60.0\" --filter-name \"FS_gt60\" \\
                       --filter-expression \"MQ < 40.0\" --filter-name \"MQ_lt40\" \\
                       --filter-expression \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum_ltm8\" \\
                       --filter-expression \"AF < 0.05\" --filter-name \"AF_05\" \\
                       --filter-expression \"AF > 0.95\" --filter-name \"AF_95\"\n";
		print OUT1 "$controlArgs{GATK} SelectVariants -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.abf.vcf -select-type SNP --restrict-alleles-to BIALLELIC -O $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf\n"; 

######STEP FOUR: Output new supercontig FASTA with ambiguity codes
		print "STEP 4: Generating IUPAC FASTA file\n";
		print OUT1 "$controlArgs{GATK} FastaAlternateReferenceMaker -R $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta -O $controlArgs{IUPAC_OUT}/$tax/$tax.iupac.fasta -V $controlArgs{GENOTYPE_OUT}/$tax/$tax.snps.biallelic.vcf --use-iupac-sample $controlArgs{GENOTYPE_OUT}/$tax/$tax.ref.fasta\n";

######STEP FIVE: Split BAM and VCF by locus then phase each locus
		print "STEP 5: Splitting BAM and VCF file by locus\n";
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
		print "STEP 6: Take break to distribute jobs\n";
		system "sbatch $controlArgs{GENOTYPE_OUT}/$tax/$tax.genotype.sh";
    }
    
    elsif ($passedArgs{runmode} == 2)
    {
    	print "STEP 7: Building phased haplotype sequences for $tax\n";
		if ($passedArgs{genotype} == 1)
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
		    if (exists $refSeqs{$tax}{$locus})
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
								if ($passedArgs{genotype} == 0)
								{
									$thisSeq = $refSeqs{$tax}{$locus};
								}
								elsif ($passedArgs{genotype} == 1)
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
								    		if ($passedArgs{genotype} == 0)
											{
								    			$phasedSeq = "$phasedSeq" . "N";
								    		}
								    		elsif ($passedArgs{genotype} == 1)
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


if ($passedArgs{runmode} == 2)
{
    print "STEP 8: Output fasta files of orthologous phased haplotype sequences\n";
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
							    $switch = 1;
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

	foreach my $locus (sort keys %locusList)
	{
    	open OUT1,'>',"$controlArgs{FASTA_OUT}/PHASED/$locus.fasta";
    	foreach my $tax (sort keys %taxaPloidy)
    	{
    		if ($passedArgs{genotype} == 0)
    		{
				if (exists $refSeqs{$tax}{$locus})
				{
				    if (exists $tax2phase{$locus}{$tax}[0])
	    			{
						for my $i (0..(scalar(@{$tax2phase{$locus}{$tax}}) - 1))
						{
					    	print OUT1 ">$tax2phase{$locus}{$tax}[$i]\n$phasedSeq{$locus}{$tax2phase{$locus}{$tax}[$i]}\n";
						}
		    		}
	    			else
				    {
			    		my $referenceHeader = "$tax" . "__REF";
						print OUT1 ">$referenceHeader\n$refSeqs{$tax}{$locus}\n";
	    			}
				}
			}
			if ($passedArgs{genotype} == 1)
    		{
				if (exists $altSeqs{$tax}{$locus})
				{
				    if (exists $tax2phase{$locus}{$tax}[0])
	    			{
						for my $i (0..(scalar(@{$tax2phase{$locus}{$tax}}) - 1))
						{
					    	print OUT1 ">$tax2phase{$locus}{$tax}[$i]\n$phasedSeq{$locus}{$tax2phase{$locus}{$tax}[$i]}\n";
						}
		    		}
	    			else
				    {
			    		my $referenceHeader = "$tax" . "__REF";
						print OUT1 ">$referenceHeader\n$altSeqs{$tax}{$locus}\n";
	    			}
				}
			}
    	}
    	close OUT1;
	}
	
########
#Added per-locus genotype fasta files to FASTA_OUT 20201109
########	
	print "STEP 9: Output fasta files of orthologous unphased genotype sequences\n";
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
	print "STEP 10: Collecting Phasing Statistics\n";
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
			if (exists $refSeqs{$tax}{$locus})
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
					$het = $nvar/(length($refSeqs{$tax}{$locus}));
				}
				$phasein{$tax}{$locus}{nvar} = $nvar;
				$phasein{$tax}{$locus}{het} = $het;
				$phasein{$tax}{$locus}{nblocks} = $nblocks;
				$phasein{$tax}{$locus}{lbl} = $lbl;
    		}
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
			if (exists $refSeqs{$tax}{$locus})
			{
				my $refLen = length($refSeqs{$tax}{$locus});
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
				$globalStats{$tax}{loclen} = $globalStats{$tax}{loclen} + length($refSeqs{$tax}{$locus});
				
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
			elsif (! exists $refSeqs{$tax}{$locus})
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