#!/bin/bash

module load gcc
module load cmake
module load boost/.1.71.0
module load mathematica/12

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#	Variable Declaration
#__________________________________________________________________________________________
#
Selection="$1"
ResultsDir=/"$2"
RunDir="$3"
Ion1="$4"
Ion2="$5"
RedThick="$6"
Fluct="$7"
CrossSec="$8"
Begin="$9"
End=${10}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#	Relevant Directories
#__________________________________________________________________________________________
#
#	CHANGE THE PROGRAM AND SHELL DIRECTORIES TO YOUR LOCATION
TrentoDir=/projects/jnorhos/pcarzon/TrentoGluonSpectrum/build_UIUC/src
ICCINGDir=/projects/jnorhos/pcarzon/ICCING
AnalysisDir=/projects/jnorhos/pcarzon/ICCING/NewAnalysis
ShellDir=/projects/jnorhos/pcarzon/ICCINGShells
DataType=$Ion1$Ion2$RedThick
DataDir=$ResultsDir/$RunDir/"$DataType"_K"$Fluct"_CS"$CrossSec"
DataFile="$DataType"_K"$Fluct"_CS"$CrossSec"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#?????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# 	USE THIS TO SUBMIT JOBS: KEYWORD = Submit
#?????????????????????????????????????????????????????????????????????????????????????????????????????????????????


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#	TRENTO Event Submission
#__________________________________________________________________________________________
#
if [[ "$Selection" == "All" || "$Selection" == "Submit" ]]
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	This if statement controls what scripts are run.
#	All = Submit Events, Combine Events, and do Calculations.
#	Comb = Combine Events and do Calculations.
#	Calc = Just do Calculations.
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
then

	mkdir $ResultsDir/$RunDir
	mkdir $DataDir
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#	Create Output directories
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	cd $ShellDir

	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#	Take general TRENTO config file and substitute run params.
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#	How sed is used here.
#	'sed' is a bash command that can substitute phrases in a file, i.e. replace all instances of NUM with "42"
#	There are many ways to use sed and I have chosen to use an easily read way:
#		- sed -e "s|ION1|$Ion1|" ...
#		- The '-e' tells sed that you are making multiple substitutions.
#		- Using the "" tells sed to interput all special characters such as '$'.
#		- The character that follows 's' is the seperator, any character will do as long as it is not in
#		  the variable being past i.e. don't use 's/' when trying to pass a file path. There must be 3
#		  of these seperators per statement.
#		- After the first seperator is the string to be replaced in the file.
#		- After the second seperator is the variable to be substituted in the file.
#	The end of the command is 'GENERAL.file > temp.file', take the general file, replace strings and write
#	as temp file.
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#	Copy TRENTO config file to data directory to give context to data.
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	for ((i=$Begin; i<=$End; i+=1));
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#	Loop to submit TRENTO events seperated among several jobs.
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		do
			mkdir $DataDir/"$i"
			sed -e "s|OUTPUT_DIR|$DataDir/"$i"|" -e "s|ION1|$Ion1|" -e "s|ION2|$Ion2|" -e "s|P_|$RedThick|" -e "s|K_|$Fluct|" -e "s|CROSS_SECTION|$CrossSec|" TRENTO_GEN.conf > TrentoParameters"$i".conf
			sed -e "s|PDIR|$TrentoDir|" -e "s|DDIR|$DataDir/"$i"|" -e "s|CONFFILE|$ShellDir/TrentoParameters"$i".conf|" -e "s|DTYPE|$DataType|" -e "s|NUM|$i|" submitTrento.sh > "$DataFile"_"$i".sh;
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			#	Generate TRENTO submission script.
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

			chmod +x "$DataFile"_"$i".sh;
			chmod 755 "$DataFile"_"$i".sh;
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			#	Change permissions of files.
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			#	'chmod' changes the permissions of given files.
			#	'+x' make file executable
			#	'755' is complicated and accomplishes the same thing as '+x'
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			TrentoID=$(sbatch -p qgp -A qgp ./"$DataFile"_"$i".sh);
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			#	qsub: submit job to UIUC Campus Cluster
			#	'-q qgp' submit job to 'qgp' queue
			#	'-l nodes=1' request access to 1 node
			#	'ppn=1' request access to 1 'cpu core'
			#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			sed -e "s|INPUTDIR|$DataDir/"$i"/|" -e "s|OUTPUTDIR|$DataDir/"$i"/|" ICCING_GEN.conf > ICCINGParameters_"$i".conf
			ICCINGID=$(sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ICCINGParameters_$i.conf",ICCINGDIR="$ICCINGDir" --dependency=afterany:${TrentoID##* } ./submitICCING.sh);
	#		sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ICCINGParameters_$i.conf" ./submitICCING.sh;

			sbatch -p qgp -A qgp --export=All,DataDir="$DataDir/$i/" --dependency=afterany:${ICCINGID##* } ./clean.sh
		done
		wait
#	mv TrentoParameters.conf ICCINGParameters_0.conf $DataDir;

fi



















#???????????????????????????????????????????????????????????????????????????????????????????????????????????????
#	USE THIS TO COMBINE FILES: KEYWORD = Comb
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#	Combine Events
#__________________________________________________________________________________________
#
if [[ "$Selection" == "Test"  ||  "$Selection" == "Comb" ]]
then
	cd $ShellDir

	mv TrentoParameters.conf ICCINGParameters_0.conf $DataDir;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="quark_counts",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="energy_eccentricities",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="baryon_eccentricities_neg",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="baryon_eccentricities_pos",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="strange_eccentricities_neg",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="strange_eccentricities_pos",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="charge_eccentricities_neg",Begin="$Begin",End="$End" ./combine2.sh;
	sbatch -p qgp -A qgp --export=All,DataDir="$DataDir",DataFile="charge_eccentricities_pos",Begin="$Begin",End="$End" ./combine2.sh;
fi


#???????????????????????????????????????????????????????????????????????????????????????????????????????????????
#	USE THIS TO RUN ANALYSIS: KEYWORD = Analyze
#???????????????????????????????????????????????????????????????????????????????????????????????????????????????

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#	Combine Events
#__________________________________________________________________________________________
#
if [[ "$Selection" == "Analyze" ]]
then
	cd $ShellDir

	mkdir $DataDir/Observables/;
	mkdir $DataDir/Observables/energy/;
	mkdir $DataDir/Observables/baryon_neg/;
	mkdir $DataDir/Observables/baryon_pos/;
	mkdir $DataDir/Observables/strange_neg/;
	mkdir $DataDir/Observables/strange_pos/;
	mkdir $DataDir/Observables/charge_neg/;
	mkdir $DataDir/Observables/charge_pos/;

	sed -e "s|INPUTFILE|$DataDir/energy_eccentricities.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/energy/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_energy.conf;
	sed -e "s|INPUTFILE|$DataDir/baryon_eccentricities_neg.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/baryon_neg/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_baryon_neg.conf;
	sed -e "s|INPUTFILE|$DataDir/baryon_eccentricities_pos.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/baryon_pos/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_baryon_pos.conf;
	sed -e "s|INPUTFILE|$DataDir/strange_eccentricities_neg.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/strange_neg/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_strange_neg.conf;
	sed -e "s|INPUTFILE|$DataDir/strange_eccentricities_pos.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/strange_pos/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_strange_pos.conf;
	sed -e "s|INPUTFILE|$DataDir/charge_eccentricities_neg.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/charge_neg/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_charge_neg.conf;
	sed -e "s|INPUTFILE|$DataDir/charge_eccentricities_pos.dat|" -e "s|OUTPUTFOLDER|$DataDir/Observables/charge_pos/|" /NewAnalysis/ANALYSIS_GEN.conf > ANALYSISParameters_charge_pos.conf;

	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_energy.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/energy/" ./submitANALYSIS.sh;
	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_baryon_neg.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/baryon_neg/" ./submitANALYSIS.sh;
	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_baryon_pos.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/baryon_pos/" ./submitANALYSIS.sh;
	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_strange_neg.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/strange_neg/" ./submitANALYSIS.sh;
	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_strange_pos.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/strange_pos/" ./submitANALYSIS.sh;
	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_charge_neg.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/charge_neg/" ./submitANALYSIS.sh;
	sbatch -p qgp -A qgp --export=All,ConfigFile="$ShellDir/ANALYSISParameters_charge_pos.conf",ANALYSISDir="$AnalysisDir",OUTPUTFOLDER="$DataDir/Observables/charge_pos/" ./submitANALYSIS.sh;
fi
