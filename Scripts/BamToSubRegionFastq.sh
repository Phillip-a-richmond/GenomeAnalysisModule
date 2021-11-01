#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 60 Gb of RAM for the whole job
#SBATCH --mem=160G

## Using 16 CPUs
#SBATCH --cpus-per-task=20

## Running for a max time of 48 hours
#SBATCH --time=48:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

##########
# Set up #
##########

NSLOTS=20
REFGENOME=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
PICARDJAR=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/opt/BamsurgeonEnvironment/share/picard-2.23.8-0/picard.jar
GENEBREAKER=/mnt/common/WASSERMAN_SOFTWARE/GeneBreaker/
cd /mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/

# Dad
# Activate BamSurgeon
BAMSURGEON=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/
source $BAMSURGEON/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BAMSURGEON/opt/BamsurgeonEnvironment

# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Father_subregion.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Father_subregion.namesorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Father_R1.fastq
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Father_R2.fastq
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR1955499_GRCh38.sorted.bam
REGION=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/chr19_stk11.bed

# Samtools extract subregion
samtools view -L $REGION $INBAM -b -o $OUTBAM 
samtools index $OUTBAM

# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED

# Convert to fastq
samtools fastq \
	-@ $NSLOTS \
	-1 $OUTFQ1 -2 $OUTFQ2 \
	-0 /dev/null \
	-s /dev/null \
	-n $OUTBAM_SORTED

#########################3

# Mom
# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Mother_subregion.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Mother_subregion.namesorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Mother_R1.fastq
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Mother_R2.fastq
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR1955435_GRCh38.sorted.bam

# Samtools extract subregion
samtools view -L $REGION $INBAM -b -o $OUTBAM
samtools index $OUTBAM

# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED

# Convert to fastq
samtools fastq \
        -@ $NSLOTS \
        -1 $OUTFQ1 -2 $OUTFQ2 \
        -0 /dev/null \
        -s /dev/null \
        -n $OUTBAM_SORTED

######################

# Proband 
# Define variables
OUTBAM=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Proband_subregion.bam
OUTBAM_SORTED=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Proband_subregion.namesorted.bam
OUTFQ1=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Proband_R1.fastq
OUTFQ2=/mnt/common/WASSERMAN_SOFTWARE/bamsurgeon/Proband_R2.fastq
INBAM=/mnt/common/OPEN_DATA/POLARIS_PROCESS/ERR2304565_GRCh38.sorted.bam

# Samtools extract subregion
samtools view -L $REGION $INBAM -b -o $OUTBAM
samtools index $OUTBAM

# Samtools sort & index
samtools sort -n $OUTBAM -@ $NSLOTS -o $OUTBAM_SORTED

# Convert to fastq
samtools fastq \
        -@ $NSLOTS \
        -1 $OUTFQ1 -2 $OUTFQ2 \
        -0 /dev/null \
        -s /dev/null \
        -n $OUTBAM_SORTED

