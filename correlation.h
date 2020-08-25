#ifndef Correlator_H
#define Correlator_H

//#undef __STRICT_ANSI__

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include <functional>

using namespace std;

class IO;

class Correlator
{
private:

  string dipole_model;
  double gamma;

  double GeVfm = 1/0.19732687;
  double R = 0.61803399;
  double C = 1.0 - R;
	void CopyCorrelator(const Correlator &e);

public:

  Correlator(string model);
  ~Correlator();

  Correlator(const Correlator &original);
  Correlator& operator=(const Correlator& original);

  double Vaccum(double r, double alpha, double m, double Qs);

  double MVModel(double r, double alpha, double m, double Qs);

  double FindMaximum(double alpha, double m, double Qs, double lower, double upper, double tolerance);

  double F(double r = 0., double alpha = 0., double m = 0., double Qs = 0.);
};
#endif
