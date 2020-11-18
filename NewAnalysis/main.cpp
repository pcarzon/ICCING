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
		inOut.OutputObservables(calculate_centbins.CalculateCummulants(), "/" + centrality_cuts[cent] + "per_centbins/Cummulants_4and6part.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_2Particle_Cummulants(), "/" + centrality_cuts[cent] + "per_centbins/Cummulants_2part.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_4and6Particle_Cummulants(), "/" + centrality_cuts[cent] + "per_centbins/Cummulants_4and6part_noRatio.dat");
		inOut.OutputObservables(calculate_centbins.Calculate_V2_V3_CummulantRatio(), "/" + centrality_cuts[cent] + "per_centbins/Cummulants_E2_E3_Ratio.dat");
	}
/*	Sorter sort_events(events);
	vector<vector<Event>> events_1per_centbins = sort_events.SortEccentricitiesIntoCentralityBins(1);
	vector<vector<Event>> events_5per_centbins = sort_events.SortEccentricitiesIntoCentralityBins(5);
	vector<vector<Event>> events_10per_centbins = sort_events.SortEccentricitiesIntoCentralityBins(10);

	Calculator calculate_1per_centbins(events_1per_centbins);
	inOut.OutputObservables(calculate_1per_centbins.CalculateCummulants(), "/1per_centbins/Cummulants_4and6part.dat");
	inOut.OutputObservables(calculate_1per_centbins.Calculate_2Particle_Cummulants(), "/1per_centbins/Cummulants_2part.dat");
	inOut.OutputObservables(calculate_1per_centbins.Calculate_4and6Particle_Cummulants(), "/1per_centbins/Cummulants_4and6part_noRatio.dat");
	inOut.OutputObservables(calculate_1per_centbins.Calculate_V2_V3_CummulantRatio(), "/1per_centbins/Cummulants_E2_E3_Ratio.dat");
//	inOut.OutputObservables(calculate_1per_centbins.CalculatePtObservables(), "/1per_centbins/Cummulants_PtObservables.dat");

	Calculator calculate_5per_centbins(events_5per_centbins);
	inOut.OutputObservables(calculate_5per_centbins.CalculateCummulants(), "/5per_centbins/Cummulants_4and6part.dat");
	inOut.OutputObservables(calculate_5per_centbins.Calculate_2Particle_Cummulants(), "/5per_centbins/Cummulants_2part.dat");
	inOut.OutputObservables(calculate_5per_centbins.Calculate_4and6Particle_Cummulants(), "/5per_centbins/Cummulants_4and6part_noRatio.dat");
	inOut.OutputObservables(calculate_5per_centbins.Calculate_V2_V3_CummulantRatio(), "/5per_centbins/Cummulants_E2_E3_Ratio.dat");
//	inOut.OutputObservables(calculate_5per_centbins.CalculatePtObservables(), "/5per_centbins/Cummulants_PtObservables.dat");

	Calculator calculate_10per_centbins(events_10per_centbins);
	inOut.OutputObservables(calculate_10per_centbins.CalculateCummulants(), "/10per_centbins/Cummulants_4and6part.dat");
	inOut.OutputObservables(calculate_10per_centbins.Calculate_2Particle_Cummulants(), "/10per_centbins/Cummulants_2part.dat");
	inOut.OutputObservables(calculate_10per_centbins.Calculate_4and6Particle_Cummulants(), "/10per_centbins/Cummulants_4and6part_noRatio.dat");
	inOut.OutputObservables(calculate_10per_centbins.Calculate_V2_V3_CummulantRatio(), "/10per_centbins/Cummulants_E2_E3_Ratio.dat");
//	inOut.OutputObservables(calculate_10per_centbins.CalculatePtObservables(), "/10per_centbins/Cummulants_PtObservables.dat");
*/
	return 0;
}
