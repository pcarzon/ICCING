#ifndef SORTER_H
#define SORTER_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################
#include "event.h"
//__________________________________________________________________________________________

using namespace std;

// compareByLength - compare by multiplicity
bool compareByEntropy(const Event &a, const Event &b)
{
    return a.entropy > b.entropy;
}

// compareByNpart -
bool compareByNpart(const Event &a, const Event &b)
{
    return a.number_of_participants > b.number_of_participants;
}

class Sorter
{
private:
  //******************************************************************************************
  //  Internal variables
  //******************************************************************************************
  vector<Event> sorted_events;
  int events_1percent_centrality_bin = 0;
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Internal Functions
//##########################################################################################
  //  Copy function for Sorter class, called by operator= and implicit copy functions
	void CopySorter(const Sorter &e);

  // GenerateCentralityCuts - sort in centrality (make file of cuttoffs)
/*  double GenerateCentralityCuts(string folder, string sort_file, vector<eccs> all_eccentricities);

  // GenerateNumberOfParticipantsCuts - sort in npart (make file of cuttoffs)
  int GenerateNumberOfParticipantsCuts(string folder, vector<eccs> all_eccentricities);
  */
//__________________________________________________________________________________________

public:

//__________________________________________________________________________________________
//##########################################################################################
//  Basic Class Functions
//##########################################################################################
  Sorter(vector<Event> all_events);  //  Class constructor, must specify path to configuration file
  ~Sorter();  //  Destructor, used to clear all data stored by class

  Sorter(const Sorter &original); //  Implicit copy function, newSorterObject(oldSorterObject)
  Sorter& operator=(const Sorter& original);  //  Defines what happens when you use = operator on class
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Sorter Specific Functions
//##########################################################################################

// sort - sorting events into centrality bins (actually does the work) (sorted_eccentricities[centrality][i])
void SortEccentricitiesIntoCentralityBins(int bin);
//__________________________________________________________________________________________
};
#endif
