#!/bin/bash
#SBATCH --job-name=__RUNID__
#SBATCH --output=__LOGFILE__
#SBATCH --mail-user=<Your Email>
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=<Your Queue>
#SBATCH --account=<Your Account>

module load picard/2.9.2
module load bwa/0.7.17
module load gatk/4.1.4.0
module load samtools/1.10
module load bamtools/2.1.1
