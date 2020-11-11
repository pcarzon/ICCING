#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>

#include "io.h"
#include "event.h"

using namespace std;


int main (int argc, char *argv[])
{
	IO inOut(argv[1]);

	vector<Event> events = inOut.ReadEvents();

	return 0;
}
