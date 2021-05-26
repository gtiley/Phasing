#!/usr/bin/perl -w
use strict;

#----------------------------------------------------------------------------------------#
#George P. Tiley and Andrew A. Crowl
#25 May 2021
#contact: george.tiley@duke.edu
#contact: andrew.crowl@duke.edu
#prepare reference sequences and reads for PATÉ
#----------------------------------------------------------------------------------------#

# Accepted AND necessary commands
my @checkArgs = ("filetype","inputFolder","outputFolder","ploidyFile");
my %passedArgs = ();
if (scalar(@ARGV) == 0)
{
die "/*--------INPUT PARAMETERS--------*/\n
--filetype STRING <fasta or fastq>
--inputFolder STRING <name of folder with input fasta or fastq files>
--outputFolder STRING <name of output folder of fasta or fastq files formatted for PATÉ>
--ploidyFile STRING <name of ploidy file>

\n/*--------EXAMPLE COMMAND--------*/\n
    perl PATE_formatInput.pl --filetype fasta --inputFolder supercontigs --outputFolder referenceSequences --ploidyFile ploidy.txt\n
    perl PATE_formatInput.pl --filetype fastq --inputFolder readsWithLongNames --outputFolder rawReads --ploidyFile ploidy.txt\n

\n/*--------FLAG OPTIONS--------*/\n
--filetype
	fasta = Supercontig output from hybpiper is expected as input
	fastq = uses the ploidy.txt file to get individual names and decompress and rename fastq files as IND.R1.fq an IND.R2.fq
	
\n/*--------NOTES--------*/\n
Meant to be a convient tool for users to go straight from hybpiper into PATÉ. Other pre-processing tools can be added upon request, but consult the example input for formatting if not using hybpiper supercontigs as reference sequences.\n";
}
elsif (scalar(@ARGV) > 0)
{
    for my $i (0..(scalar(@ARGV) - 1))
    {
		if ($ARGV[$i] eq "--filetype")
		{
	    	$passedArgs{filetype} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "--inputFolder")
		{
	    	$passedArgs{inputFolder} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "--outputFolder")
		{
        	$passedArgs{outputFolder} = $ARGV[$i+1];
        	system "mkdir $passedArgs{outputFolder}";
		}
		if ($ARGV[$i] eq "--ploidyFile")
		{
	    	$passedArgs{ploidyFile} = $ARGV[$i+1];
		}
	}
	foreach my $arg (@checkArgs)
	{
		if (! exists $passedArgs{$arg})
		{
	    	die "/*--------MISSING PARAMETER--------*/\nMissing command line argument: $arg\n\n";
		}
	}
}

if ($passedArgs{filetype} eq "fasta")
{
	my %locusList = ();
	my %seqs = ();
	my %taxList = ();
	open FH1,'<',"$passedArgs{ploidyFile}";
	while (<FH1>)
	{
		if (/^(\S+)\s+.+/)
		{
			my $tax = $1;
			if (! exists $taxList{$tax})
			{
				$taxList{$tax} = 1;
			}
		}
	}
	close FH1;

	my @ff1 = glob("$passedArgs{inputFolder}/*.fasta");
	foreach my $ff (@ff1)
	{
    	if ($ff =~ m/$passedArgs{inputFolder}\/(\S+)\.fasta/)
    	{
        	my $locus = $1;
        	my $tax = "";
	        if (! exists $locusList{$locus})
    	    {
        	    $locusList{$locus} = 1;
        	}
       		open FH1,'<',"$ff";
	        while (<FH1>)
    	    {                                                                                                         
            	if (/^>(\S+)/)
            	{
                	$tax = $1;
					foreach my $ind (sort keys %taxList)
					{
						if (index($tax,$ind) >= 0)
						{
							$tax = $ind;
						}
					}
                	$seqs{$locus}{$tax} = "";
            	}
	            elsif (/(\S+)/)
    	        {
        	        my $seq = $1;
            	    $seqs{$locus}{$tax} = $seqs{$locus}{$tax} . $seq;
            	}
        	}
	        close FH1;
    	}
	}

	foreach my $locus (sort keys %locusList)
	{
    	open OUT1,'>',"$passedArgs{outputFolder}/$locus.fasta";
	    foreach my $tax (sort keys %taxList)
    	{
        	if (exists $seqs{$locus}{$tax})
	        {
    	        print OUT1 ">$tax\n$seqs{$locus}{$tax}\n";
        	}
    	}
    	close OUT1;
	}
}


if ($passedArgs{filetype} eq "fastq")
{
	my %taxList = ();
	open FH1,'<',"$passedArgs{ploidyFile}";
	while (<FH1>)
	{
		if (/^(\S+)\s+.+/)
		{
			my $tax = $1;
			if (! exists $taxList{$tax})
			{
				$taxList{$tax} = 1;
			}
		}
	}
	close FH1;
	
	foreach my $tax (sort keys %taxList)
	{
		my @fqFiles = glob("$passedArgs{inputFolder}/*$tax*.*");
		foreach my $ff (@fqFiles)
		{
    		if ($ff =~ m/$passedArgs{inputFolder}\/\S*R1\S*\.(\S+)/)
    		{
	        	my $extension = $1;
	        	if ($extension =~ m/gz/)
	        	{
    	    		system "cp $ff $passedArgs{outputFolder}/$tax.R1.fq.gz";
    	    	}
    	    	elsif ($extension =~ m/fq/ || $extension =~ m/fastq/)
	        	{
    	    		system "cp $ff $passedArgs{outputFolder}/$tax.R1.fq";
    	    	}
    	    }
    	    elsif ($ff =~ m/$passedArgs{inputFolder}\/\S*R2\S*\.(\S+)/)
    		{
	        	my $extension = $1;
	        	if ($extension =~ m/gz/)
	        	{
    	    		system "cp $ff $passedArgs{outputFolder}/$tax.R2.fq.gz";
    	    	}
    	    	elsif ($extension =~ m/fq/ || $extension =~ m/fastq/)
	        	{
    	    		system "cp $ff $passedArgs{outputFolder}/$tax.R2.fq";
    	    	}
    	    }
    	}
	}
}
exit;