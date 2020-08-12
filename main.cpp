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
	ofstream quark_output;
	if (inOut.GetTest() == "QuarkRatio")
	{	quark_output.open(inOut.GetOutputDir() + "quark_ratio_test.dat");	}

	while (!inOut.LastEvent())
	{
		start = clock();
		testEvent = inOut.ReadEvent(initializedEvent);
		cout << "started processing event" << endl;
			while (!testEvent.IsEventDone())
			{
				Sample testSample;
				Quarks testQuarks;

				testSample = testEvent.SampleEnergy();

				if (testSample.q_s == -100)
				{
					break;
				}
				testQuarks = machine.SplitSample(testSample);
				if (inOut.GetTest() == "QuarkRatio")
				{	quark_output << testSample.q_s << " " << testQuarks.GetCharge()[0] << endl;	}

				testEvent.UpdateDensity(testQuarks);
				eventcount++;
			}
		cout << "# times through event loop " << eventcount << endl;
		eventcount = 0;
		inOut.WriteEvent(testEvent);

		if (inOut.GetTest() == "QuarkRatio")
		{	quark_output.close();	}

		duration = (clock() - start)/(double)CLOCKS_PER_SEC;
		cout << "Event processing time: " << duration/60 << " min" << endl;
	}

	return 0;
}
