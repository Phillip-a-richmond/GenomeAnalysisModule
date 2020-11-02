#!/bin/bash

#SBATCH --partition=defq

#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#SBATCH --array=0-8%1

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
WORKING_DIR=/mnt/common/OPEN_DATA/POLARIS_PROCESS/
mkdir -p $WORKING_DIR
RAW_DIR=/mnt/common/OPEN_DATA/POLARIS_RAW/
NSLOTS=$SLURM_CPUS_PER_TASK

# Define sample id from set of fastq files, based on the job array index
Files=(${RAW_DIR}*_1.fastq.gz)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SampleR1Fastq=${array[-1]}
IFS='_' read -a array2 <<< "${SampleR1Fastq}"
SAMPLE=${array2[0]}

echo $SAMPLE
echo "${SAMPLE}_1.fastq.gz"
echo "${SAMPLE}_2.fastq.gz"

FASTQR1=$RAW_DIR${SAMPLE}_1.fastq.gz
FASTQR2=$RAW_DIR${SAMPLE}_2.fastq.gz

ls $FASTQR1
ls $FASTQR2

# STEP 1: MAP AGAINST GRCH37 #

# Genome Information
#### GRCh37 ####
BWA_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa
# Set SAMPLE_ID to have GRCh37
SAMPLE_ID=${SAMPLE}_GRCh37

##################
# Initialize Log #
##################
# Timing from here to end of script
FullRunStart=`date +%s`

LOGFILE=${SAMPLE_ID}_logfile.csv
rm ${SAMPLE_ID}_logfile.csv
touch ${SAMPLE_ID}_logfile.csv
# This logfile will be written to in csv format.
# Columns for log file:
echo "SampleID,Operation,Runtime">> $LOGFILE

#############
# Map reads #
#############
cd $WORKING_DIR

Start=`date +%s`

bwa mem $BWA_INDEX \
        -t $NSLOTS \
        -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" \
        -M \
        $FASTQR1 \
        $FASTQR2 \
        > $WORKING_DIR$SAMPLE_ID.sam

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,BWAmem,$runtime" >> $LOGFILE
echo "BWA mem ran in $runtime"


############
# Samtools #
############
#convert to binary and index
Start=`date +%s`

samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
        | samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,Samtools,$runtime" >> $LOGFILE
echo "Samtools ran in $runtime"

rm $WORKING_DIR$SAMPLE_ID'.sam'


# STEP 1: MAP AGAINST GRCH38 #

# Genome Information
#### GRCh38 ####
BWA_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GRCh38-lite.fa
# Set SAMPLE_ID to have GRCh38
SAMPLE_ID=${SAMPLE}_GRCh38

##################
# Initialize Log #
##################
# Timing from here to end of script
FullRunStart=`date +%s`

LOGFILE=${SAMPLE_ID}_logfile.csv
rm ${SAMPLE_ID}_logfile.csv
touch ${SAMPLE_ID}_logfile.csv
# This logfile will be written to in csv format.
# Columns for log file:
echo "SampleID,Operation,Runtime">> $LOGFILE

#############
# Map reads #
#############
cd $WORKING_DIR

Start=`date +%s`

bwa mem $BWA_INDEX \
        -t $NSLOTS \
        -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" \
        -M \
        $FASTQR1 \
        $FASTQR2 \
        > $WORKING_DIR$SAMPLE_ID.sam

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,BWAmem,$runtime" >> $LOGFILE
echo "BWA mem ran in $runtime"


############
# Samtools #
############
#convert to binary and index
Start=`date +%s`

samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
        | samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,Samtools,$runtime" >> $LOGFILE
echo "Samtools ran in $runtime"

rm $WORKING_DIR$SAMPLE_ID'.sam'

exit
