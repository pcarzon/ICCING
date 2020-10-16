#!/bin/bash

#source ../../../../ect/profile.d/module.sh

#module load gcc
#module load cmake
#module load boost/1.71.0
#module load mathematica/12

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#	Variable Declaration
#__________________________________________________________________________________________
#

#SBATCH --mem-per-cpu=4450

ConfigFile=$1

cd /projects/jnorhos/pcarzon/ICCING;

chmod +x submitICCING.sh;
chmod 755 submitICCING.sh;
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	Change permissions of files.
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#	'chmod' changes the permissions of given files.
#	'+x' make file executable
#	'755' is complicated and accomplishes the same thing as '+x'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
make iccing;
sbatch -p qgp --export=All,ConfigFile="run_parameters.conf" ./submitICCING.sh;
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	qsub: submit job to UIUC Campus Cluster
#	'-q qgp' submit job to 'qgp' queue
#	'-l nodes=1' request access to 1 node
#	'ppn=1' request access to 1 'cpu core'
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
