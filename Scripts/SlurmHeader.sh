#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=yourEmailHere@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 60 Gb of RAM for the whole job
#SBATCH --mem=60G

## Using 16 CPUs
#SBATCH --cpus-per-task=16

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

