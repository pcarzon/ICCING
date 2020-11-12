#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <string>
#include <string.h>
#include <vector>
#include <iostream>


using namespace std;

struct cent{
// ways to define centrality
//cls centrality class?
	double cls, multiplicity, number_of_participants;
	};


struct Event
	{
		// si is entropy
		// M multiplicity
		// pT average pT
		// MID multiplicity of identified particles
		// psi for flow usually
		// phi for eccs
		// ND uses no decays
		int event_num = 0;
		int number_of_participants = 0;
		double entropy = 0, impact_parameter = 0;
		double multiplicity = 0, average_pt = 0, radius = 0;
		double eccentricity[7], phi[7];
		double flow_harmonics[7], psi[7], flow_harmonics_no_decays[7], psi_no_decays[7];

	};


#endif
