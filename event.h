#ifndef Event_H
#define Event_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################
#include "global.h"
#include "ecc.h"
//__________________________________________________________________________________________

using namespace std;

//__________________________________________________________________________________________
//##########################################################################################
//  Class Declarations
//##########################################################################################
class IO;
//__________________________________________________________________________________________
extern default_random_engine get_random_number;

class Event
{
private:
  //__________________________________________________________________________________________
  //##########################################################################################
  //  Event Input Parameters
  //##########################################################################################
  double gluon_rad;
  double quark_rad;
  double lambda_;
  double grid_max;
  double grid_step;
  double tau_0;
  double e_thresh;
  double kappa_;  //  Used for Qs grid
  int grid_points;
  string test_;
  uniform_int_distribution<int> get_grid_point;
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Density Grids
  //##########################################################################################
  vector<vector<double>> initial_energy;  //  Input Energy density
  vector<vector<double>> t_a; //  Target Input Energy density
  vector<vector<double>> t_b; //  Projectile Input Energy density
  vector<vector<vector<double>>> density; //  ICCING densities: 0(gluon), 1(baryon), 2(em_charge), 3(strange), 4(charm)
  vector<vector<int>> gluon_dist;  //  Sample from initial_energy for ICCING algorithm
  vector<vector<double>> quark_dist; //  Projectile Input Energy density
  vector<vector<int>> valued_points;
  vector<vector<vector<double>>> eccentricities;
  Eccentricity ecc;
  int number_gluon = 0;
  int number_up = 0;
  int number_down = 0;
  int number_strange = 0;
  int number_charm = 0;
  int x_center;
  int y_center;

  double total_initial_energy = 0;  //  Records initial total of initial_energy
  double total_initial_entropy = 0;
  double total_energy = 0;  //  Keeps track of the total of initial_energy as ICCING runs
  double total_entropy = 0;
  Sample out_sample;
  double seed;  //  Records the random seed used in event
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Internal Functions
  //##########################################################################################
	//  Copy function for Event class, called by operator= and implicit copy functions
  void CopyEvent(const Event &e);

  //  Select energy of gluon
  Sample GetGlue();

  //  Gets intigration bounds for density grid manipulations
  vector<int> GetIntegrationBounds(int size, double raduis, double xcenter, double ycenter);

  //  Subtracts energy from initial_energy and adds it to density[0]
  void UpdateEnergy(double ratio);
  //__________________________________________________________________________________________

public:

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Basic Class Functions
  //##########################################################################################
  Event();  // Class Constructor
  ~Event(); //  Class Destructor

  Event(const Event &original); //  Implicit copy function, newIOObject(oldIOObject)
  Event& operator=(const Event& original);  //  Defines what happens when you use = operator on class
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Event Specific Functions
  //##########################################################################################
  //  Sample Initial Energy for ICCING algorithm
  Sample SampleEnergy();

  //  Propogates Results of Splitter
  bool UpdateDensity(Quarks quark_density);

  //  Calculate Eccentricities from density grids
  void CalculateEccentricities();

  //  Clears Event variables as a cautionary measure
  void CleanEvent();

  //  Checks Event totals and returns true when initial_total is below a threshold
  bool IsEventDone();
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Friend Classes
  //##########################################################################################
  //  IO class can access private variables of Event class
  friend class IO;
  //__________________________________________________________________________________________
};
#endif
