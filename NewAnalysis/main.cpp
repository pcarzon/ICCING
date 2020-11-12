#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>

#include "io.h"
#include "event.h"
#include "sorter.h"

using namespace std;


int main (int argc, char *argv[])
{
	IO inOut(argv[1]);

	vector<Event> events = inOut.ReadEvents();

	Sorter Centrality_10percent_bins(events);
	Centrality_10percent_bins.SortEccentricitiesIntoCentralityBins(10);

	return 0;
}
