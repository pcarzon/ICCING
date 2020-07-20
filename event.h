#ifndef Event_H
#define Event_H

//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <ofstream>
#include <string>
#include <cmath>
#include <vector>
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
  double e_thresh;
  int grid_points;
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

  double total_initial_energy;  //  Records initial total of initial_energy
  double total_energy;  //  Keeps track of the total of initial_energy as ICCING runs
  double seed;  //  Records the random seed used in event
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Internal Functions
  //##########################################################################################
	//  Copy function for Event class, called by operator= and implicit copy functions
  void CopyEvent(const Event &e);

  //  See: RollGlue in ICCING_v0_1_8.nb
  double RollGlue(double e_tot);

  //  See: RollGlue in ICCING_v0_1_8.nb
  double GetQs(int x_center, int y_center);

  //  See: RollGlue in ICCING_v0_1_8.nb
  double GetEtot(int x_center, int y_center);
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
  void UpdateDensity(Quarks quark_density);

  //  Clears Event variables as a cautionary measure
  void CleanEvent();
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
