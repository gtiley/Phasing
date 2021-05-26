#input and output paths
REF = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/referenceSequences
GENOTYPE_OUT = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/genotypeOutput2
PHASE_OUT = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/phasedOutput2
FASTA_OUT = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/fastaOutput3
IUPAC_OUT = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/iupacOutput2
SUMMARYSTATS_OUT = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/summaryStatsOutput
FQ = /blue/burleigh/gtiley/phasing/20200412/Dryopteris/Dryopteris_RawReads
PLOIDY = /blue/burleigh/gtiley/phasing/20200412/Dryopteris_hybpiper/ploidy.txt

#paths to software or just the binaries if already in your path
BWA = bwa 
PICARD = picard
GATK = gatk
SAMTOOLS = samtools
BAMTOOLS = bamtools
BGZIP = bgzip
TABIX = tabix
HPOPG = /blue/burleigh/gtiley/phasing/H-PoPG_2/H-PoPG/H-PoPGv0.2.0.jar
SCHEDULER = sbatch

#gatk filter options
GATK_FILTER_EXPRESSION = "QD < 2.0" "QD_lt2"
GATK_FILTER_EXPRESSION = "FS > 60.0" "FS_gt60"
GATK_FILTER_EXPRESSION = "MQ < 40.0" "MQ_lt40"
GATK_FILTER_EXPRESSION = "ReadPosRankSum < -8.0" "ReadPosRankSum_ltm8"
GATK_FILTER_EXPRESSION = "AF < 0.05" "AF_05"
GATK_FILTER_EXPRESSION = "AF > 0.95" "AF_95"