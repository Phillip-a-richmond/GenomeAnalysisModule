#!/bin/bash

########################################################################
###
### Copyright (c) 2020, Wasserman lab
###
### FILE        install.sh
###
### DESCRIPTION This is the install.sh script for the project.
###             We use conda in order to install the application.
###
### Initial version @ Godfrain Jacques KOUNKOU
### Modified @ Phillip Richmond
########################################################################

# set -e


##########################################
# DESC This function will find all FlexTyper dependencies and build FlexTyper 
# ARGS This function doesnt require any arguments
# RSLT The side effect is the newly created FlexTyper directory with the FlexTyper binary
##########################################
function Bamsurgeon() 
{
	DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
	OPT_DIR=${DIR}/opt
	mkdir -p ${OPT_DIR}
	pushd ${OPT_DIR}
	MINI_CONDA_INSTALL_DIR=$OPT_DIR/miniconda3

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINI_CONDA_INSTALL_DIR 

	source ${MINI_CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
	conda env create --prefix ${OPT_DIR}/BamsurgeonEnvironment -f $DIR/Bamsurgeon_Downstream_CondaEnv.yml 
	cd $DIR
	conda activate opt/BamsurgeonEnvironment
}


Bamsurgeon
