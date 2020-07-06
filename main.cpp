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

int main (int argc, char *argv[])
{
	IO inputOutputObject(argv[1]);
/*"/projects/jnorhos/pcarzon/ICCING/testInput/run_parameters.conf"*/
	Event testEvent;
	testEvent.ReadInitialEnergy(inputOutputObject.ReadEvent());
	inputOutputObject.WriteEvent(testEvent);
	cout << "Hello World!" << endl;

	return 0;
}
