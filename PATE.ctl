#input and output paths
REF = YOUR_PATH/referenceSequences
GENOTYPE_OUT = YOUR_PATH/genotypeOutput
PHASE_OUT = YOUR_PATH/phasedOutput
FASTA_OUT = YOUR_PATH/fastaOutput
IUPAC_OUT = YOUR_PATH/iupacOutput
SUMMARYSTATS_OUT = YOUR_PATH/summaryStatsOutput
FQ = YOUR_PATH/rawReads
PLOIDY = YOUR_PATH/ploidy.txt

#paths to software or just the binaries if already in your path
BWA = bwa 
PICARD = picard
GATK = gatk
SAMTOOLS = samtools
BAMTOOLS = bamtools
BGZIP = bgzip
TABIX = tabix
HPOPG = YOUR_PATH/H-PoPG_2/H-PoPG/H-PoPGv0.2.0.jar
SCHEDULER = sbatch

#gatk filter options
GATK_FILTER_EXPRESSION = "QD < 2.0" "QD_lt2"
GATK_FILTER_EXPRESSION = "FS > 60.0" "FS_gt60"
GATK_FILTER_EXPRESSION = "MQ < 40.0" "MQ_lt40"
GATK_FILTER_EXPRESSION = "ReadPosRankSum < -8.0" "ReadPosRankSum_ltm8"
GATK_FILTER_EXPRESSION = "AF < 0.05" "AF_lt05"
GATK_FILTER_EXPRESSION = "AF > 0.95" "AF_gt95"
GATK_FILTER_EXPRESSION = "DP < 10" "DP_lt10"