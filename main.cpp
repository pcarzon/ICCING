#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>

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
	int eventcount = 0;
	clock_t start;
	double duration;

	initializedEvent = inOut.InitializeEvent();
	inOut.InitializeEOS();

	while (!inOut.LastEvent())
	{
		start = clock();
		testEvent = inOut.ReadEvent(initializedEvent);
		cout << "started processing event" << endl;
			while (!testEvent.IsEventDone())
			{

		testSample = testEvent.SampleEnergy();

//		testQuarks = machine.SplitSample(testSample);

//		testEvent.UpdateDensity(testQuarks);
		eventcount++;
			}
			cout << "# times through event loop " << eventcount << endl;
			eventcount = 0;
		inOut.WriteEvent(testEvent);

		duration = (clock() - start)/(double)CLOCKS_PER_SEC;
		cout << "Event processing time: " << duration/60 << " min" << endl;
	}

	return 0;
}
