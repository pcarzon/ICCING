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
	initializedEvent = inOut.InitializeEvent();
	inOut.InitializeEOS();

//	while (!inOut.LastEvent())
	//{
		testEvent = inOut.ReadEvent(initializedEvent);
		Sample testSample = testEvent.SampleEnergy();
		machine.SplitSample(testSample);
		inOut.WriteEvent(testEvent);
	//}
	cout << "Hello World!" << endl;

	return 0;
}
