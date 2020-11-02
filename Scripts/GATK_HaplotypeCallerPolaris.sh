#!/bin/bash

#SBATCH --partition=defq

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#SBATCH --array=1-8%1

##########
# Set up #
##########

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Get GATK in the path
# unload_bcchr
source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
# load_cvmfs
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
# load gatk
module load gatk/3.7 

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR

# Set working space
WORKING_DIR=/mnt/common/OPEN_DATA/POLARIS_VCF/
mkdir -p $WORKING_DIR
BAM_DIR=/mnt/common/OPEN_DATA/POLARIS_PROCESS/
NSLOTS=$SLURM_CPUS_PER_TASK

# Get the files we're going to use
Files=(${BAM_DIR}*.sorted.bam)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SAMPLE_BAM=${array[-1]}
IFS='_' read -a array2 <<< "${SAMPLE_BAM}"
SAMPLE=${array2[0]}
GENOME=${array2[1]::-11}

echo $SAMPLE
echo $GENOME
ls $BAM_DIR/$SAMPLE_BAM

if [ $GENOME = "GRCh37" ];
then
        #### GRCh37 ####
        echo "GRCh37 genome"
        BWA_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa
        SAMPLE_ID=${SAMPLE}_GRCh37
	DBSNP=/mnt/common/DATABASES/REFERENCES/GRCh37/DBSNP/All_20180423.vcf
fi

if [ $GENOME = "GRCh38" ];
then

        #### GRCh38 #### 
        echo "GRCh38 genome"
        BWA_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
        SAMPLE_ID=${SAMPLE}_GRCh38
	DBSNP=/mnt/common/DATABASES/REFERENCES/GRCh38/DBSNP/All_20180423.vcf
fi

echo $SAMPLE_ID
echo $BWA_INDEX

## Running GATK
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller \
	-nct $NSLOTS \
	-R $BWA_INDEX \
	-I $BAM_DIR/$SAMPLE_BAM \
	--pcr_indel_model NONE \
	-o $WORKING_DIR/${SAMPLE_ID}_GATK_HCv3.7.vcf

# And GVCF mode
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller \
	-nct $NSLOTS \
        -R $BWA_INDEX \
	--emitRefConfidence GVCF \
        -I $BAM_DIR/$SAMPLE_BAM \
        --pcr_indel_model NONE \
        -o $WORKING_DIR/${SAMPLE_ID}_GATK_HCv3.7.g.vcf


