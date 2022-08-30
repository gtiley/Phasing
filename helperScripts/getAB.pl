#!/usr/bin/perl -w
use strict;

my $vcfFile = $ARGV[0];
my $outFileRoot = $ARGV[1];
my %ABdata = ();
my %nvariants = ();
my %chromosomes  = ();
my %positions = ();
my @taxa = ();
my $ntax = 0;
open FH1,'<',"$vcfFile";
my $skippingHeader = 1;
while (<FH1>)
{
	my $line =$_;
	chomp $line;
	if ($line =~ m/^#CHROM/)
	{
	    my @temp = ();
	    @temp = split(/\s+/,$line);
	    for my $i (9..(scalar(@temp)-1))
	    {
		my @refPath = split(/\//,$temp[$i]);
		$taxa[$i] = $refPath[(scalar(@refPath)-2)];
		$nvariants{$taxa[$i]} = 0;
		$ntax++;
	    }
	    $skippingHeader = 0;
	}
	if ($skippingHeader == 0)
	{
		my @temp = ();
		@temp = split(/\s+/,$line);
		if ($temp[6] eq "PASS")
		{
			for my $i (9..(scalar(@temp)-1))
			{
			    my $filterValues = $temp[$i];
			    my @filterValuesVector = ();
			    @filterValuesVector = split(/\:/,$filterValues);
			    #print "@filterValuesVector\n";
			    if ($filterValuesVector[1] =~ m/(\d+)\,(\d+)/)
			    {
				my $refAllele = $1;
				my $altAllele = $2;
				if ((($refAllele + $altAllele) >= 10) && ($refAllele > 1) && ($altAllele > 1))
				{
				    my $abValue = $altAllele/($refAllele + $altAllele);
				    push @{$ABdata{$taxa[$i]}}, $abValue;
				    push @{$chromosomes{$taxa[$i]}}, $temp[0];
				    push @{$positions{$taxa[$i]}}, $temp[1];
				    $nvariants{$taxa[$i]} = $nvariants{$taxa[$i]} + 1;
				}
			    }
			}
		}
	}
}
close FH1;

for my $i (9..($ntax + 8))
{
    open OUT1,'>',"$outFileRoot/$taxa[$i]/$taxa[$i].ab";
    print OUT1 "Chromosome\tPosition\tAB\n";
    for my $j (0..($nvariants{$taxa[$i]} - 1))
    {
	print OUT1 "$chromosomes{$taxa[$i]}[$j]\t$positions{$taxa[$i]}[$j]\t$ABdata{$taxa[$i]}[$j]\n";
    }
    close OUT1;
}
exit;
