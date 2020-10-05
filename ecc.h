#ifndef Eccentricity_H
#define Eccentricity_H

#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;
// Do these systems PbPb, XeXe, OO(woods saxon), pPb (nucleon width = 0.3)
// 3 million events each
class Eccentricity
{
private:

  double x_center_of_mass = 0;
  double y_center_of_mass = 0;
  vector<vector<double>> sparse_density;

	void CopyEccentricity(const Eccentricity &e);
  void CleanEccentricity();

  //vector<double> Eccentricities(vector<vector<double>> grid, double grid_step);

  vector<double> StandardCalculation(string density_type, int m, int n);

  vector<double> NewCalculation(string density_type, int m, int n);

public:

  Eccentricity();
  ~Eccentricity();

  Eccentricity(const Eccentricity &original);
  Eccentricity& operator=(const Eccentricity& original);

  vector<vector<vector<double>>> CalculateEccentricities(int grid_max, double grid_step, vector<vector<vector<double>>> density);
};
#endif
