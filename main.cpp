#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "ecc.h"
#include "eos.h"
#include "event.h"
#include "io.h"
#include "splitting.h"

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
/*"/projects/jnorhos/pcarzon/ICCING/testInput/run_parameters.conf"*/
	Splitter machine = inOut.InitializeSplitter();

	Event testEvent, initializedEvent;

	Sample testSample;
	Quarks testQuarks;

	initializedEvent = inOut.InitializeEvent();
	inOut.InitializeEOS();

//	while (!inOut.LastEvent())
	//{
		testEvent = inOut.ReadEvent(initializedEvent);

		//	while (!testEvent.IsEventDone())
			//{

		testSample = testEvent.SampleEnergy();

		testQuarks = machine.SplitSample(testSample);

		testEvent.UpdateDensity(testQuarks);
		//	}

		inOut.WriteEvent(testEvent);
	//}

	return 0;
}
