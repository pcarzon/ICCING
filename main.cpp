//__________________________________________________________________________________________
//##########################################################################################
//	ICCING Code Version 1.0
//	Designed by Matthew Sievert
// 	Implimented and Evolved by Patrick Carzon
// 	Written at Rutgers University and University of Illinois at Urbana-Champaign
//
//	For questions and comments please email:
//	Patrick Carzon @ pcarzon2@illinois.edu
//##########################################################################################
//__________________________________________________________________________________________


//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################
#include "ecc.h"
#include "event.h"
#include "io.h"
#include "splitting.h"
#include "correlation.h"
//__________________________________________________________________________________________

using namespace std;

//__________________________________________________________________________________________
//##########################################################################################
//  Random Number Generators
//##########################################################################################
default_random_engine get_random_number;
//__________________________________________________________________________________________

// 	When running a code from the command line and passing values to it, the information
//	is read in as argc, which is the number of strings passed to the program with the
//	first one being the name of the program by default, and argv, which is an array of
//	the strings passed as arguments to the program
int main (int argc, char *argv[])
{
	//	Declare and initialize the input/output object
	//	This gets passed the path to the config file and reads in operating parameters
	IO inOut(argv[1]);

	//******************************************************************************************
  //  Declare relevent variables and objects
  //******************************************************************************************
	Event currentEvent, initializedEvent;
	Splitter machine;

	//******************************************************************************************
  //  Initialize objects
  //******************************************************************************************
	//	Set parameters relevent to the splitter class
	//	Also read in the flavor chemistry file
	machine = inOut.InitializeSplitter();

	//	Set parameters relevent to the event class
	//	Also initialize density profiles for sampling gluons and distributing quarks
	initializedEvent = inOut.InitializeEvent();

	//	Read in and process the equation of state file
	inOut.InitializeEOS();

	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// 	This is a test that records the mass of quark pairs and their relevent Qs
	ofstream quark_output;
	if (inOut.GetTest() == "QuarkRatio"){	quark_output.open(inOut.GetOutputDir() + "quark_ratio_test.dat");	}
	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//******************************************************************************************
  //  Main Event Loop
  //******************************************************************************************
	while (!inOut.LastEvent())
	{
		//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		//	Read next event using initializedEvent as base
		currentEvent = inOut.ReadEvent(initializedEvent);
		//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			//******************************************************************************************
	  	//  Event Loop, Process event until initial energy density is empty
	  	//******************************************************************************************
			while (!currentEvent.IsEventDone())
			{
				bool successful_density_update;
				//	Declare Sample and Quarks for individual event processing
				Sample currentSample;
				Quarks currentQuarks;

				//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				//	Get an energy sample from event
				currentSample = currentEvent.SampleEnergy();
				//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

				// If initial energy density is less than e_thresh, continue event loop and find new point
				if (currentSample.q_s == -100){	continue;	}

				//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				//	Generate Quarks from event energy sample
				currentQuarks = machine.SplitSample(currentSample);
				//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

				//	If there was not enough energy to create 2 quarks of given flavor mass, sample event again
				if (currentQuarks.GetEnergyFraction() == -1)	{	continue;	}

				//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// 	This is a test that records the mass of quark pairs and their relevent Qs
				if (inOut.GetTest() == "QuarkRatio") {	quark_output << currentSample.q_s << " " << currentQuarks.GetCharge()[0] << endl;	}
				//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				//	Update Density grids with sampled quarks
				successful_density_update = currentEvent.UpdateDensity(currentQuarks);
				//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

				// 	If quarks are out of bounds of grid, sample event again
				if (!successful_density_update) { cout << "Grid is too small. Quarks out of bounds." << endl;	exit(0); }

			}

		//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		//	Calculate Eccentricities of event
		currentEvent.CalculateEccentricities();

		//	Write event data to files
		inOut.WriteEvent(currentEvent);
		//	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	}

	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// 	This is a test that records the mass of quark pairs and their relevent Qs
	if (inOut.GetTest() == "QuarkRatio"){	quark_output.close();	}
	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	return 0;
}
