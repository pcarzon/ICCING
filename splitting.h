#ifndef Splitter_H
#define Splitter_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "global.h"
#include "functions.h"

using namespace std;

class IO;
extern default_random_engine get_random_number;
//extern vector<SplineSet> CubicSpline(vector<double> &x, vector<double> &y);
//extern SplineSet SplineInterval(double value);
//extern double InterpolateValue(SplineSet poly, double value);

class Splitter
{
private:

	vector<vector<SplineSet>> flavor_chemistry;  //  flavor_chemistry[Qs(GeV)][quark_prob] quark_prob: 0 = Qs, 1 = up, 2 = down ...
	string dipole_model;  //  Future Update
	double alpha_s;
	double alpha_min;
	double r_max;
	double e_thresh;
	double lambda_;

	//  See: RollGlue in ICCING_v0_1_8.nb
  double RollGlue(double e_tot);

	Charge RollFlavor(double Qs);

	void CopySplitter(const Splitter &e);

public:

  Splitter();
  ~Splitter();

  Splitter(const Splitter &original);
  Splitter& operator=(const Splitter& original);

  Quarks SplitSample(Sample sampled_energy);  //  See: SplittingSample in ICCING_v0_1_8.nb

	friend class IO;
};
#endif
