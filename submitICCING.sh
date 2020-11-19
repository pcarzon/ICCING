#!/bin/bash

#source /ect/profile.d/modules.sh

#SBATCH --mem-per-cpu=8450

module load gcc/7.2.0
module load cmake/3.12.0
module load boost/1.71.0

cd $ICCINGDIR

./iccing $ConfigFile
