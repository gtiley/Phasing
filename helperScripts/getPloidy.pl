#!/usr/bin/perl -w
use strict;

my $phasingRoot = $ARGV[0];
my $ploidyFile = $ARGV[1];
my $estimatedPloidyPath = $ARGV[2];

my %ploidies = ();

open FH1,'<',"$phasingRoot/$ploidyFile";
while(<FH1>)
{
	if (/^(\S+)\s+(\d+)\s+.*/)
	{
		my $tax = $1;
		my $startPloidy = $2;
		if ($startPloidy != 2)
		{
			print "Warning - starting ploidy for estimation was not 2 ($tax = $startPloidy) and results may be flawed\n";
		}
		$ploidies{$tax} = $startPloidy;
		open FH2,'<',"$estimatedPloidyPath/$tax/$tax.modelSelection.txt";
		while(<FH2>)
		{
			if (/Ploidy\:\s+(\d+)/)
			{
				my $estimatedPloidy = $1;
				$ploidies{$tax} = $estimatedPloidy;
			}
		}
		close FH2;
	}
}
close FH1;

open OUT1,'>',"$phasingRoot/$ploidyFile.estimated";
foreach my $tax (keys %ploidies)
{
	print OUT1 "$tax\t$ploidies{$tax}\n"
}
close OUT1;