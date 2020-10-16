#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <string>
#include <string.h>
#include <vector>
#include <iostream>

using namespace std;

struct pars{

	int m,n,p;
	int fm,fn,fp;

	int rat;

	int swit;
	};
	
struct cent{
// ways to define centrality
	double cls,M,npart;
	};

struct eccs
	{
		int npart;
		double si,b;
		double v[7],M,pT,R,MID;
		double psi[7],vND[7],psiND[7];
		double ec[7],phi[7];
	};

struct SCsub
	{
		double vnvm,m4,v4,v6,m6;
	};


struct serr
	{
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

