#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 60 Gb of RAM for the whole job
#SBATCH --mem=160G

## Using 16 CPUs
#SBATCH --cpus-per-task=40

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

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment


# Variables are first defined (no spaces allowed!)
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/ProblemSet4/
# Make a problem set 4 directory
mkdir -p $WORKING_DIR
# Change to it
cd $WORKING_DIR

# Arrange your data
cp /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/Files/Sample1*fastq $WORKING_DIR

# Change this if necessary to match your sample ID
SAMPLE=Sample1
FASTQR1=${SAMPLE}_R1.fastq
FASTQR2=${SAMPLE}_R2.fastq

GENOME_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa

####################
# Map with BWA mem #
####################

bwa mem $GENOME_INDEX \
	-R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:illumina" \
	$FASTQR1 \
	$FASTQR2 \
	> ${SAMPLE}.sam

######################################
# Convert, sort, index with samtools #
######################################

samtools view \
	-b ${SAMPLE}.sam  \
	-o ${SAMPLE}.bam

samtools sort \
	${SAMPLE}.bam \
	-o ${SAMPLE}.sorted.bam

samtools index \
	${SAMPLE}.sorted.bam

##################################
# Call variants with Singularity #
##################################
# Load singularity
module load singularity
BIN_VERSION="1.0.0"

# Pull latest version, if you already have it, this will be skipped
singularity pull docker://google/deepvariant:"${BIN_VERSION}"

# I'm still using the same SAMPLE defined above
SAMPLE_BAM=$SAMPLE.sorted.bam

# Then, I set up the reference genome location
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/
FASTA_FILE=GRCh38-lite.fa

# Last, we'll make a special output directory called ${SAMPLE_ID}_Deepvariant
OUTPUT_DIR=$WORKING_DIR/${SAMPLE}_DeepVariant
mkdir -p $OUTPUT_DIR

# This is a complicated command, which will work if you set the above correctly.
singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
	-B "${WORKING_DIR}":"/bamdir" \
	-B "${FASTA_DIR}":"/genomedir" \
	-B "${WORKING_DIR}":"/output" \
	docker://google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/genomedir/$FASTA_FILE" \
  --reads="/bamdir/$SAMPLE_BAM" \
  --output_vcf="/output/${SAMPLE}_DeepVariant.vcf.gz" \
  --num_shards=$SLURM_CPUS_PER_TASK 


