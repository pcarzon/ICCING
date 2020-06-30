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
	IO inputOutputObject;

	Event testEvent;

	testEvent.SampleEnergy(1, "Hello!");
	cout << "Hello World!" << endl;

	return 0;
}
