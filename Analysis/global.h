#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <string>
#include <string.h>
#include <vector>
#include <iostream>

using namespace std;

struct pars{

	// p momentum bin?
	int m,n,p;
	int fm,fn,fp;

	int rat;

	int swit;
	};

struct cent{
// ways to define centrality
//cls centrality class?
	double cls, multiplicity, number_of_participants;
	};

struct eccs
	{
		// si is entropy
		// M multiplicity
		// pT average pT
		// MID multiplicity of identified particles
		// psi for flow usually
		// phi for eccs
		// ND uses no decays
		int number_of_participants;
		double entropy, impact_parameter;
		double flow_harmonics[7], multiplicity, average_pt, radius;
		double psi[7], flow_harmonics_no_decays[7], psi_no_decays[7];
		double eccentricity[7], phi[7];
	};

struct SCsub
	{
		// symetric cummulants without normalization
		double vnvm,m4,v4,v6,m6;
	};


struct serr
	{
		//
		double cen,ans,ans2,err,err2;
	};

inline string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) elems.push_back(item);
    }
    return elems;
}


inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

#endif
