#!/bin/bash

# Here I have a few variables, which I define as VARNAME="<someVariable,like a path>"
WorkshopDir=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/

ls $WorkshopDir
ls /mnt/scratch/Public/TRAINING/GenomeAnalysisModule/

WorkingDir=/mnt/scratch/Public/TRAINING/GenomeAnalysisModule/StudentSpaces/Sherlock/

cd $WorkingDir

# There are also environmental variables
# $PWD
echo $PWD

echo "I am in the directory: $PWD"

Name=Sherlock
echo "Hi, my name is $Name!"

# You can make variables within variables by using ${}
MyWorkshopDir=/scratch/tr-precisionhealth-1/Workshops/StudentSpaces/${Name}

