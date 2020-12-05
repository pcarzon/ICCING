#ifndef GLOBAL_H_
#define GLOBAL_H_

//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <string>
#include <string.h>
#include <vector>
#include <iostream>
//__________________________________________________________________________________________

using namespace std;

//	Data structure for events
struct Event
	{
		int event_num = 0;
		int number_of_participants = 0;
		double entropy = 0, impact_parameter = 0;
		double multiplicity = 0, average_pt = 0, radius = 0;
		double eccentricity[7], phi[7];
		double flow_harmonics[7], psi[7], flow_harmonics_no_decays[7], psi_no_decays[7];
	};

#endif
