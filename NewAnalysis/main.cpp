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

	vector<int> centrality_cuts = {1, 5, 10};

	for (int cent = 0; cent < centrality_cuts.size(); cent++)
	{
		Sorter sort_events(events);
		vector<vector<Event>> events_centbins = sort_events.SortEccentricitiesIntoCentralityBins(centrality_cuts[cent]);

		Calculator calculate_centbins(events_centbins);
		inOut.MakeDirectory(centrality_cuts[cent]);
		inOut.OutputObservables(calculate_centbins.CalculateCummulants(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_4and6part.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_2Particle_Cummulants(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_2part.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_4and6Particle_Cummulants(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_4and6part_noRatio.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_V2_V3_CummulantRatio(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_E2_E3_Ratio.dat");
	}

	return 0;
}
