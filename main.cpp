#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "ecc.h"
#include "eos.h"
#include "event.h"
#include "io.h"
#include "functions.h"
#include "splitting.h"

using namespace std;
//class Event;
//class IO;
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
	//Event testEvent;
	Splitter machine = inOut.InitializeSplitter();

	Event testEvent, initializedEvent;
	initializedEvent = inOut.InitializeEvent();

//	while (!inOut.LastEvent())
	//{

		testEvent = inOut.ReadEvent(initializedEvent);
		machine.SplitSample(testEvent.SampleEnergy());
		inOut.WriteEvent(testEvent);
	//	testEvent.CleanEvent();
	//	(&testEvent)->~Event();
	//	new (&testEvent) Event();
	//}
	cout << "Hello World!" << endl;

	return 0;
}
