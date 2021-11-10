#!/bin/bash

#SBATCH --partition=training_q

## Change to be your email address
#SBATCH --mail-user=yourEmailHere@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 80 Gb of RAM for the whole job
#SBATCH --mem=80G

## Using 10 CPUs
#SBATCH --cpus-per-task=10

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

## Get the tools we need, from a conda environment 
source  /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda  activate  GenomeAnalysis


