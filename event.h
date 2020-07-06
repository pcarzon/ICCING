#ifndef Event_H
#define Event_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "io.h"

using namespace std;

class IO;
class Eccentricity;

class Event
{
private:

  vector<vector<double>> initial_energy;
  vector<vector<double>> t_a;
  vector<vector<double>> t_b;
  vector<vector<vector<double>>> density;
  vector<vector<double>> sampled_energy;

  double total_initial_energy;
  double total_energy;
  double seed;

	void CopyEvent(const Event &e);

public:

  Event();
  ~Event();

  Event(const Event &original);
  Event& operator=(const Event& original);

  void ReadInitialEnergy(vector<vector<double>> initEnergy);
  void SampleEnergy(int location, double radius);
  void UpdateDensity(int densityType, vector<vector<double>> inputDensity);

  vector<vector<double>> GetInitialEnergy() { return initial_energy; }
  vector<vector<double>> GetTa() { return t_a; }
  vector<vector<double>> GetTb() { return t_b; }

  friend vector<vector<double>> IO::ReadEvent();
};
#endif
