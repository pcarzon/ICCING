#ifndef Event_H
#define Event_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "global.h"
using namespace std;

class IO;
class Eccentricity;

class Event
{
private:

  double grid_max;
  double grid_step;
  int grid_points;

  vector<vector<double>> initial_energy;
  vector<vector<double>> t_a;
  vector<vector<double>> t_b;
  vector<vector<vector<double>>> density;
  vector<vector<double>> sampled_energy;

  double total_initial_energy;
  double total_energy;
  double seed;

	void CopyEvent(const Event &e);

  double RollGlue(double e_tot);  //  See: RollGlue in ICCING_v0_1_8.nb
public:

  Event();
  ~Event();

  Event(const Event &original);
  Event& operator=(const Event& original);

  void ReadInitialEnergy(vector<vector<double>> initEnergy);
  Sample SampleEnergy(); //  See: First 2 commands in While in DistributeCharge in ICCING_v0_1_8.nb
    //  Calls RollGlue
  void UpdateDensity(Quarks quark_density);

  vector<vector<double>> GetInitialEnergy() { return initial_energy; }
  vector<vector<double>> GetTa() { return t_a; }
  vector<vector<double>> GetTb() { return t_b; }



  void CleanEvent();

  friend class IO;
};
#endif
