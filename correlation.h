#ifndef Correlator_H
#define Correlator_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include <functional>
//__________________________________________________________________________________________

using namespace std;

class IO;

class Correlator
{
private:
  //__________________________________________________________________________________________
  //##########################################################################################
  //  Correlator Parameters
  //##########################################################################################
  // Flag for type of correlation function used in quark position and momentum selection
  string dipole_model;
  double lambda_bym;
  double alpha_s;

  // Function data type which stores a reference to correlation function being used
  function<double(double, double, double, double)> corr;

  double GeVfm = 1/0.19732687;  // Conversion factor between GeV and fm
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Internal Functions
  //##########################################################################################
	//  Copy function for Correlator class, called by operator= and implicit copy functions
	void CopyCorrelator(const Correlator &e);

  //  Correlation function in GBWModel
  double GBWModel(double r, double alpha, double m, double Qs);

  //  Correlation function in MV Model
  double MVModel(double r, double alpha, double m, double Qs);
  //__________________________________________________________________________________________

public:
  //__________________________________________________________________________________________
  //##########################################################################################
  //  Basic Class Functions
  //##########################################################################################
  Correlator(string model, double lambda, double alpha_s);  // Class Constructor
  ~Correlator();  //  Class Destructor

  Correlator(const Correlator &original); //  Implicit copy function, newIOObject(oldIOObject)
  Correlator& operator=(const Correlator& original);  //  Defines what happens when you use = operator on class
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Internal Functions
  //##########################################################################################
  //  Find Maximum Probability given quark flavor, random alpha, and Qs
  //  Golden Section Search
  double FindMaximum(double alpha, double m, double Qs, double lower, double upper, double tolerance);

  // Given quark pair data, determine distance and momentum fraction
  double F(double r = 0., double alpha = 0., double m = 0., double Qs = 0.);
  //__________________________________________________________________________________________
};
#endif
