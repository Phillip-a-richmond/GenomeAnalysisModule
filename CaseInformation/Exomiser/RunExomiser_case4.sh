#!/bin/bash
  
#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 8 Gb of RAM for the whole job
#SBATCH --mem=8G

## Using 1 CPUs
#SBATCH --cpus-per-task=1

## Running for a max time of 20 minutes
#SBATCH --time=00:20:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


#Case 4
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/mnt/common/WASSERMAN_SOFTWARE/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/Exomiser/Case4_Exomiser.yml

cd /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/Exomiser/
mkdir -p /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/Exomiser/Case4/
cd /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/CaseInformation/Exomiser/Case4/

java -Xmx4g -jar -Djava.io.tmpdir=$PWD $EXOMISER_DIR/$EXOMISER \
         --analysis $CONFIG_YML \
        --spring.config.location=$EXOMISER_DIR
