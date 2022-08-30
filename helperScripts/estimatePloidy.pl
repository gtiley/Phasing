#!/usr/bin/perl -w

$scheduler = $ARGV[0];
$template = $ARGV[1];
$ploidypath = $ARGV[2];
$helperpath = $ARGV[3];
$tax = $ARGV[4];

$scheduler =~ s/scheduler\_//;

open OUT1,'>',"$ploidypath/$tax/$tax.mixturemodels.sh";
open FH1,'<',"$template";
while(<FH1>)
{
my $line = $_;
chomp $line;
$line =~ s/__RUNID__/$tax.mixturemodels.id/;
$line =~ s/__LOGFILE__/$tax.mixturemodels.log/;
print OUT1 "$line\n";
}
close FH1;
print OUT1 "$ploidypath/$tax\n";
open OUT2,'>',"$ploidypath/$tax/getPloidy.R";
print OUT2 "setwd(\"$ploidypath/$tax\");\n";
print OUT2 "source(\"$helperpath/Ks_plots/ploidy.test.R\")\;\nsource(\"$helperpath/Ks_plots/fitMixEM.R\")\;\nsource(\"$helperpath/Ks_plots/plotComponentExpectations.R\")\;\n";
print OUT2 "dat <- read.table(\"$ploidypath/$tax/$tax.ab\",header=TRUE,sep=\"\\t\")\;\n";
print OUT2 "ploidy.test(dat\$AB,maxPloidy=6,model=4,nstarts=100,outPrefix=\"$tax\")\;\n";
print OUT2 "quit();\n";
close OUT2;
print OUT1 "R CMD BATCH $ploidypath/$tax/getPloidy.R\n";
close OUT1;
system "$scheduler $ploidypath/$tax/$tax.mixturemodels.sh";