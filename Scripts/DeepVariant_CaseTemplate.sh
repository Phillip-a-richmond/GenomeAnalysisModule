#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
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

# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Gaku/CaseAnalysis/
OUTPUT_DIR=$WORKING_DIR
#export SINGULARITY_TMPDIR=$WORKING_DIR

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=GRCh38
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/
FASTA_FILE=GRCh38-lite.fa

SEQ_TYPE=WGS
BAM_DIR=$WORKING_DIR
Case_ID=Case3
FAMILY_ID=$Case_ID
PROBAND_ID=${Case_ID}_proband
MOTHER_ID=${Case_ID}_mother
FATHER_ID=${Case_ID}_father

PROBAND_BAM=${PROBAND_ID}.sorted.bam
FATHER_BAM=${FATHER_ID}.sorted.bam
MOTHER_BAM=${MOTHER_ID}.sorted.bam

PROBAND_VCF=${PROBAND_ID}.vcf.gz
FATHER_VCF=${FATHER_ID}.vcf.gz
MOTHER_VCF=${MOTHER_ID}.vcf.gz

PROBAND_GVCF=${PROBAND_ID}.gvcf.gz
FATHER_GVCF=${FATHER_ID}.gvcf.gz
MOTHER_GVCF=${MOTHER_ID}.gvcf.gz


# Proband 
singularity exec -e -c \
	-B "${BAM_DIR}":"/bamdir" \
	-B "${FASTA_DIR}":"/genomedir" \
	-B "${OUTPUT_DIR}":"/output" \
       -W $OUTPUT_DIR \
       $DeepVariantSIF \
  /opt/deepvariant/bin/run_deepvariant \
  --intermediate_results_dir="/output/intermediate_results_dir" \
  --model_type=WGS \
  --ref="/genomedir/$FASTA_FILE" \
  --reads="/bamdir/$PROBAND_BAM" \
  --output_vcf="/output/$PROBAND_VCF" \
  --output_gvcf="/output/$PROBAND_GVCF" \
  --num_shards=$NSLOTS 


# Mother
singularity exec -e -c \
       -B "${BAM_DIR}":"/bamdir" \
       -B "${FASTA_DIR}":"/genomedir" \
       -B "${OUTPUT_DIR}":"/output" \
       -W $OUTPUT_DIR \
       $DeepVariantSIF \
 /opt/deepvariant/bin/run_deepvariant \
  --intermediate_results_dir="/output/intermediate_results_dir" \
 --model_type=WGS \
 --ref="/genomedir/$FASTA_FILE" \
 --reads="/bamdir/$MOTHER_BAM" \
 --output_vcf="/output/$MOTHER_VCF" \
 --output_gvcf="/output/$MOTHER_GVCF" \
 --num_shards=$NSLOTS 

# Father
singularity exec -e -c \
       -B "${BAM_DIR}":"/bamdir" \
       -B "${FASTA_DIR}":"/genomedir" \
       -B "${OUTPUT_DIR}":"/output" \
       -W $OUTPUT_DIR \
       $DeepVariantSIF \
 /opt/deepvariant/bin/run_deepvariant \
  --intermediate_results_dir="/output/intermediate_results_dir" \
 --model_type=WGS \
 --ref="/genomedir/$FASTA_FILE" \
 --reads="/bamdir/$FATHER_BAM" \
 --output_vcf="/output/$FATHER_VCF" \
 --output_gvcf="/output/$FATHER_GVCF" \
 --num_shards=$NSLOTS 

#GLNexus
/mnt/common/Precision/GLNexus/glnexus_cli -c DeepVariant${SEQ_TYPE} \
        --threads $NSLOTS \
        *gvcf.gz \
        > ${FAMILY_ID}.glnexus.merged.bcf

bcftools view ${FAMILY_ID}.glnexus.merged.bcf | bgzip -c > ${FAMILY_ID}.merged.vcf.gz


