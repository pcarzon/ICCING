#ifndef Eccentricity_H
#define Eccentricity_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
//__________________________________________________________________________________________

using namespace std;

//##########################################################################################
//  Explaination of calculation and functions provided on group wiki at
//  https://wiki.illinois.edu/wiki/display/NUCLEARTHEORY/Eccentricity+Calculation
//##########################################################################################
class Eccentricity
{
private:

//__________________________________________________________________________________________
//##########################################################################################
//  Internal Parameters
//##########################################################################################
  // center of mass for initial condition, used in calculations
  double x_center_of_mass = 0;
  double y_center_of_mass = 0;

  //  Sparse density array, includes all densities, sparse array makes calculation faster
  vector<vector<double>> sparse_density;
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Internal Functions
//##########################################################################################
  //  Copy function for Eccentricity class, called by operator= and implicit copy functions
  void CopyEccentricity(const Eccentricity &e);

  //  Calculate eccentricity
  vector<double> StandardCalculation(string density_type, int m, int n);

  //  Calculate eccentricities, seperating positive and negative density values
  vector<double> NewCalculation(string density_type, int m, int n);
//__________________________________________________________________________________________

public:

//__________________________________________________________________________________________
//##########################################################################################
//  Basic Class Functions
//##########################################################################################
  Eccentricity();  //  Class constructor
  ~Eccentricity();  //  Destructor, used to clear all data stored by class

  Eccentricity(const Eccentricity &original); //  Implicit copy function, newEccentricityObject(oldEccentricityObject)
  Eccentricity& operator=(const Eccentricity& original);  //  Defines what happens when you use = operator on class
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Eccentricity Specific Functions
//##########################################################################################
  //  Calculate All eccentricities for given event
  vector<vector<vector<double>>> CalculateEccentricities(int grid_max, double grid_step, vector<vector<vector<double>>> density);

  //  Clean class
  void CleanEccentricity();
//__________________________________________________________________________________________
};
#endif
