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
//__________________________________________________________________________________________

using namespace std;

//__________________________________________________________________________________________
//##########################################################################################
//  Class Declarations
//##########################################################################################
class IO;
class Eccentricity;
//__________________________________________________________________________________________
extern default_random_engine get_random_number;

class Event
{
private:
  //__________________________________________________________________________________________
  //##########################################################################################
  //  Event Input Parameters
  //##########################################################################################
  double kappa_; //  Used for Qs grid
  double gluon_rad;
  double quark_rad;
  double lambda_;
  double grid_max;
  double grid_step;
  double tau_0;
  double e_thresh;
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
  vector<double> eccentricities;
  int number_up = 0;
  int number_down = 0;
  int number_strange = 0;
  int number_charm = 0;
  int x_center;
  int y_center;

  double total_initial_energy;  //  Records initial total of initial_energy
  double total_energy;  //  Keeps track of the total of initial_energy as ICCING runs
  Sample out_sample;
  double seed;  //  Records the random seed used in event
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Internal Functions
  //##########################################################################################
	//  Copy function for Event class, called by operator= and implicit copy functions
  void CopyEvent(const Event &e);

  //  See: RollGlue in ICCING_v0_1_8.nb
  Sample GetGlue();

  //  See: RollGlue in ICCING_v0_1_8.nb
  vector<int> GetIntegrationBounds(int size, double raduis);

  //  Copy function for Event class, called by operator= and implicit copy functions
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
  Sample SampleEnergy(); //  See: First 2 commands in While in DistributeCharge in ICCING_v0_1_8.nb
    //  Calls RollGlue

  //  Propogates Results of Splitter
  bool UpdateDensity(Quarks quark_density);

  void CalculateEccentricities();

  //  Clears Event variables as a cautionary measure
  void CleanEvent();
  //__________________________________________________________________________________________

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
