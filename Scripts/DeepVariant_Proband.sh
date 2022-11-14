#!/bin/bash
#SBATCH --partition=defq
#SBATCH --mail-user=YourEmailAddress@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=80G
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

##########
# Set up #
##########

# Load singularity
module load singularity

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
source  /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda  activate  GenomeAnalysis

# Where is the DeepVariant SIF, and set singularity cache
unset $PYTHONPATH
export SINGULARITY_CACHEDIR=$PWD
DeepVariantSIF=/mnt/common/Precision/DeepVariant/deepvariant_1.2.0.sif

# Change to your workshop directory
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/
cd $WORKING_DIR 

# Here, I'm setting up the sample I'm going to run it on
SAMPLE_ID=Proband
SAMPLE_BAM=$SAMPLE_ID.sorted.bam

# Then, I set up the reference genome location
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/
FASTA_FILE=GRCh38-lite.fa

# Last, we'll make a special output directory called ${SAMPLE_ID}_Deepvariant
OUTPUT_DIR=$WORKING_DIR/${SAMPLE_ID}_DeepVariant
mkdir -p $OUTPUT_DIR

# This is a complicated command, which will work if you set the above correctly.
singularity exec -c -e \
	-B "${WORKING_DIR}":"/bamdir" \
	-B "${FASTA_DIR}":"/genomedir" \
	-B "${WORKING_DIR}":"/output" \
	-W $WORKING_DIR \
	$DeepVariantSIF \
  /opt/deepvariant/bin/run_deepvariant \
    --intermediate_results_dir="/output/intermediate_results_dir" \
  --model_type=WGS \
  --ref="/genomedir/$FASTA_FILE" \
  --reads="/bamdir/$SAMPLE_BAM" \
  --output_vcf="/output/${SAMPLE_ID}_DeepVariant.vcf.gz" \
  --regions "19:1,200,000-1,300,000" \
  --num_shards=$SLURM_CPUS_PER_TASK 


