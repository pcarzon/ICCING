#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>

#include "io.h"
#include "event.h"
#include "sorter.h"
#include "calculator.h"

using namespace std;


int main (int argc, char *argv[])
{
	IO inOut(argv[1]);

	vector<Event> events = inOut.ReadEvents();

	Sorter sort_events(events);
	vector<vector<Event>> events_1per_centbins = sort_events.SortEccentricitiesIntoCentralityBins(1);
	vector<vector<Event>> events_5per_centbins = sort_events.SortEccentricitiesIntoCentralityBins(5);
	vector<vector<Event>> events_10per_centbins = sort_events.SortEccentricitiesIntoCentralityBins(10);

	Calculator calculate_10per_centbins(events_10per_centbins);
	return 0;
}
