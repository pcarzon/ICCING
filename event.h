#ifndef Event_H
#define Event_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

class Event
{
private:

  vector<vector<double>> initial_energy;
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
};
#endif
