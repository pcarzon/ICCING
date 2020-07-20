#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "ecc.h"
#include "eos.h"
#include "event.h"
#include "io.h"
#include "probdist.h"
#include "splitting.h"

using namespace std;
//class Event;
//class IO;
int main (int argc, char *argv[])
{
	IO inOut(argv[1]);
/*"/projects/jnorhos/pcarzon/ICCING/testInput/run_parameters.conf"*/
	//Event testEvent;
	Splitter machine = inOut.InitializeSplitter();

	Event testEvent, initializedEvent;
	initializedEvent = inOut.InitializeEvent();

	while (!inOut.LastEvent())
	{

		testEvent = inOut.ReadEvent(initializedEvent);
		testEvent.Sample();
		inOut.WriteEvent(testEvent);
	//	testEvent.CleanEvent();
	//	(&testEvent)->~Event();
	//	new (&testEvent) Event();
	}
	cout << "Hello World!" << endl;

	return 0;
}
