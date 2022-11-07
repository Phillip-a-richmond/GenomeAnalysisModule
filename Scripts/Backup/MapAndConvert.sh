#!/bin/bash

#SBATCH --partition=defq

#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

##########
# Set up #
##########

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR

# Set working space
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/
mkdir -p $WORKING_DIR
RAW_DIR=/mnt/common/OPEN_DATA/POLARIS_RAW/
NSLOTS=$SLURM_CPUS_PER_TASK

SAMPLE=Case10
echo $SAMPLE
echo "${SAMPLE}_1.fastq.gz"
echo "${SAMPLE}_2.fastq.gz"

FASTQR1=$RAW_DIR${SAMPLE}_1.fastq.gz
FASTQR2=$RAW_DIR${SAMPLE}_2.fastq.gz

ls $FASTQR1
ls $FASTQR2

# STEP 1: MAP AGAINST GRCH38 #

# Genome Information
#### GRCh38 ####
BWA_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
# Set SAMPLE_ID to have GRCh38
SAMPLE_ID=${SAMPLE}_GRCh38

#############
# Map reads #
#############
cd $WORKING_DIR

bwa mem /mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa \
        -t 40 \
        -R "@RG\tID:Case10\tSM:Case10\tPL:illumina" \
        $FASTQR1 \
        $FASTQR2 \
        > $WORKING_DIR$SAMPLE_ID.sam

############
# Samtools #
############
#convert to binary and index
Start=`date +%s`

samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
        | samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

