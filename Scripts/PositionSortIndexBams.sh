#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 60 Gb of RAM for the whole job
#SBATCH --mem=350G

## Using 16 CPUs
#SBATCH --cpus-per-task=80

## Running for a max time of 48 hours
#SBATCH --time=48:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#SBATCH --array=0-8%9

##########
# Set up #
##########

NSLOTS=78
REFGENOME=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa

# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment


cd /mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/

# Define sample id from set of fastq files, based on the job array index
BAM_DIR=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
Files=(${BAM_DIR}/Case*proband.bam)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SampleBam=${array[-1]}
IFS='_' read -a array2 <<< "${SampleBam}"
SAMPLE=${array2[0]}


# Samtools sort & index
samtools sort $SampleBam -@ $NSLOTS -o ${SAMPLE}_position.sorted.bam 
samtools index ${SAMPLE}_position.sorted.bam 


