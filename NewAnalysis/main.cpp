//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Analysis Header Files
//##########################################################################################
#include "io.h"
#include "event.h"
#include "sorter.h"
#include "calculator.h"
//__________________________________________________________________________________________

using namespace std;

int main (int argc, char *argv[])
{
	IO inOut(argv[1]);

	// List of events read from input file
	vector<Event> events = inOut.ReadEvents();

	//	List of centrality cuts, Analysis is done for each value below
	vector<int> centrality_cuts = {1, 5, 10};

	//	Loop through centrality cuts
	for (int cent = 0; cent < centrality_cuts.size(); cent++)
	{
		//	Sort events with respect to centrality bin width
		Sorter sort_events(events);
		vector<vector<Event>> events_centbins = sort_events.SortEccentricitiesIntoCentralityBins(centrality_cuts[cent]);

		//	Run calculations on sorted and binned events
		Calculator calculate_centbins(events_centbins);
		inOut.OutputObservables(calculate_centbins.CalculateCummulants(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_4and6part.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_2Particle_Cummulants(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_2part.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_4and6Particle_Cummulants(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_4and6part_noRatio.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_V2_V3_CummulantRatio(), to_string(centrality_cuts[cent]) + "per_centbins/Cummulants_E2_E3_Ratio.dat");
	}

	return 0;
}
