#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 80 Gb of RAM for the whole job
#SBATCH --mem=160G

## Using 10 CPUs
#SBATCH --cpus-per-task=20

## Running for a max time of 72 hours
#SBATCH --time=72:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

##########
# Set up #
##########

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
source  /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda  activate  GenomeAnalysis

NSLOTS=$SLURM_CPUS_PER_TASK

# Variables are first defined (no spaces allowed!)
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Gaku/CaseAnalysis/
# Make your case analysis directory
mkdir -p $WORKING_DIR
# Change to it
cd $WORKING_DIR

# Case ID
Case_ID=Case3
# Case directory location
CASE_DIRECTORY=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/CaseFiles/${Case_ID}/
# Setting genome index variable
GENOME_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa

###########
# Proband #
###########

# Proband Data
SAMPLE=${Case_ID}_proband
FASTQR1=$CASE_DIRECTORY/${SAMPLE}_R1.fastq.gz
FASTQR2=$CASE_DIRECTORY/${SAMPLE}_R2.fastq.gz


####################
# Map with BWA mem #
####################

bwa mem $GENOME_INDEX \
	-R "@RG\tID:proband\tSM:proband\tPL:illumina" \
	-t $NSLOTS \
	$FASTQR1 \
	$FASTQR2 \
	> ${SAMPLE}.sam

######################################
# Convert, sort, index with samtools #
######################################

samtools view \
	-@ $NSLOTS \
	-b ${SAMPLE}.sam  \
	-o ${SAMPLE}.bam

rm ${SAMPLE}.sam

samtools sort \
	-@ $NSLOTS \
	${SAMPLE}.bam \
	-o ${SAMPLE}.sorted.bam

samtools index \
	${SAMPLE}.sorted.bam

# Clean up and delete SAM and unsorted bam
rm ${SAMPLE}.sam 
rm ${SAMPLE}.bam

###############################################################################

##########
# Mother #
##########

# Mother Data
SAMPLE=${Case_ID}_mother
FASTQR1=$CASE_DIRECTORY/${SAMPLE}_R1.fastq.gz
FASTQR2=$CASE_DIRECTORY/${SAMPLE}_R2.fastq.gz


####################
# Map with BWA mem #
####################

bwa mem $GENOME_INDEX \
	-R "@RG\tID:mother\tSM:mother\tPL:illumina" \
	-t $NSLOTS \
	$FASTQR1 \
	$FASTQR2 \
	> ${SAMPLE}.sam

######################################
# Convert, sort, index with samtools #
######################################

samtools view \
	-@ $NSLOTS \
	-b ${SAMPLE}.sam  \
	-o ${SAMPLE}.bam

rm ${SAMPLE}.sam

samtools sort \
	-@ $NSLOTS \
	${SAMPLE}.bam \
	-o ${SAMPLE}.sorted.bam

samtools index \
	${SAMPLE}.sorted.bam

# Clean up and delete SAM and unsorted bam
rm ${SAMPLE}.sam 
rm ${SAMPLE}.bam

###############################################################################

##########
# Father #
##########

# Father Data
SAMPLE=${Case_ID}_father
FASTQR1=$CASE_DIRECTORY/${SAMPLE}_R1.fastq.gz
FASTQR2=$CASE_DIRECTORY/${SAMPLE}_R2.fastq.gz


####################
# Map with BWA mem #
####################

bwa mem $GENOME_INDEX \
	-R "@RG\tID:father\tSM:father\tPL:illumina" \
	-t $NSLOTS \
	$FASTQR1 \
	$FASTQR2 \
	> ${SAMPLE}.sam

######################################
# Convert, sort, index with samtools #
######################################

samtools view \
	-@ $NSLOTS \
	-b ${SAMPLE}.sam  \
	-o ${SAMPLE}.bam
rm ${SAMPLE}.sam

samtools sort \
	-@ $NSLOTS \
	${SAMPLE}.bam \
	-o ${SAMPLE}.sorted.bam

samtools index \
	${SAMPLE}.sorted.bam

# Clean up and delete SAM and unsorted bam
rm ${SAMPLE}.sam 
rm ${SAMPLE}.bam

###############################################################################

