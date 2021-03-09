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
#include "greenfunction.h"
//__________________________________________________________________________________________

using namespace std;

//__________________________________________________________________________________________
//##########################################################################################
//  Random Number Generators
//##########################################################################################
default_random_engine get_random_number;
//__________________________________________________________________________________________

int main (int argc, char *argv[])
{
	IO inOut(argv[1]);

	//******************************************************************************************
  //  Declare relevent variables and objects
  //******************************************************************************************
	Event testEvent, initializedEvent;
	Splitter machine;

	clock_t start;
	double duration;

	//******************************************************************************************
  //  Initialize objects
  //******************************************************************************************
	machine = inOut.InitializeSplitter();

	initializedEvent = inOut.InitializeEvent();
	inOut.InitializeEOS();

	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ofstream quark_output;
	if (inOut.GetTest() == "QuarkRatio"){	quark_output.open(inOut.GetOutputDir() + "quark_ratio_test.dat");	}
	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//******************************************************************************************
  //  Main Event Loop
  //******************************************************************************************
	while (!inOut.LastEvent())
	{
		start = clock();	//	Start Clock for timing event

		//	Read next event using initializedEvent as base
		testEvent = inOut.ReadEvent(initializedEvent);

/*		if (inOut.GetTest() == "GreensFunction")
		{
			testEvent = initializedEvent;
			cout << "initializedEvent" << endl;

		}*/
			//******************************************************************************************
	  	//  Event Loop, Process event until initial energy density is empty
	  	//******************************************************************************************
			while (!testEvent.IsEventDone())
			{
				//	Declare Sample and Quarks for individual event processing
				Sample testSample;
				Quarks testQuarks;

				//	Get an energy sample from event
				testSample = testEvent.SampleEnergy();
				// If initial energy density is empty, end event loop and start new event
				if (testSample.q_s == -100){	continue;	}

				//	Generate Quarks from event energy sample
				testQuarks = machine.SplitSample(testSample);

				//	If there was not enough energy to create 2 quarks of given flavor mass, sample event again
				if (testQuarks.GetEnergyFraction() == -1)	{	continue;	}

				//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if (inOut.GetTest() == "QuarkRatio") {	quark_output << testSample.q_s << " " << testQuarks.GetCharge()[0] << endl;	}
				//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				//	Update Density grids with sampled quarks, if quarks are out of bounds of grid, sample event again
				if (!testEvent.UpdateDensity(testQuarks)) { continue; }

			}

		//	Calculate Eccentricities of event
		testEvent.CalculateEccentricities();

		//	Write event data to files
		inOut.WriteEvent(testEvent);

		//	Clean event and print time taken to process
		testEvent.CleanEvent();
		duration = (clock() - start)/(double)CLOCKS_PER_SEC;
		cout << "Event processing time: " << duration/60 << " min" << endl;
	}

	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (inOut.GetTest() == "QuarkRatio"){	quark_output.close();	}
	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	return 0;
}
