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

int main (void)
{
	IO inputOutputObject("/projects/jnorhos/pcarzon/ICCING/testInput/000.dat", "/projects/jnorhos/pcarzon/ICCING/testOutput/000.dat", 0, 0);

	Event testEvent;

	testEvent.SampleEnergy(1, 1.0);
	cout << "Hello World!" << endl;

	return 0;
}
