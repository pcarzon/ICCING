#include "correlation.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Correlator
//##########################################################################################
Correlator::Correlator(string model, double lambda, double alphas)
{
  dipole_model = model;
  lambda_bym = lambda;
  alpha_s = alphas;
  //  Statement that ties chosen correlation function to general corr function variable
  //  Done this way so that the correlation function only has to be determined once
  if (dipole_model == "MV")
  {
    //  If using the MV Model, attach the MVModel function to corr
    corr = bind(&Correlator::MVModel, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4);
  }
  else
  {
    //  If using no particular model, attach the GBWModel correlation function to corr
    corr = bind(&Correlator::GBWModel, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4);
  }

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
Correlator::~Correlator()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Correlator::Correlator(const Correlator &original)
{
  CopyCorrelator(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Correlator Copy Function
//##########################################################################################
void Correlator::CopyCorrelator(const Correlator &e)
{
//  corr = e.corr;
  dipole_model = e.dipole_model;
  lambda_bym = e.lambda_bym;
  alpha_s = e.alpha_s;

  if (dipole_model == "MV")
  {
    //  If using the MV Model, attach the MVModel function to corr
    corr = bind(&Correlator::MVModel, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4);
  }
  else
  {
    //  If using no particular model, attach the GBWModel correlation function to corr
    corr = bind(&Correlator::GBWModel, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4);
  }
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Correlator& Correlator::operator= (const Correlator& original)
{
	CopyCorrelator(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Correlation function in GBWModel
//##########################################################################################
double Correlator::GBWModel(double r, double alpha, double m, double Qs)
{
  double term1 = alpha_s*r*(pow(m, 2)/(4*M_PI))*pow(GeVfm, 2);
  double term2 = (1 - exp(-0.25*(pow(alpha, 2) + pow(1 - alpha, 2))*pow(GeVfm*r*Qs, 2)));
  double term3 = (pow(alpha, 2) + pow(1 - alpha, 2))*pow(cyl_bessel_k(1, GeVfm*m*r), 2) + pow(cyl_bessel_k(0, GeVfm*m*r), 2);
  return term1*term2*term3;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Correlation function in MV Model
//##########################################################################################
double Correlator::MVModel(double r, double alpha, double m, double Qs)
{
  double term1 = alpha_s*r*(pow(m, 2)/(4*M_PI))*pow(GeVfm, 2);
  double term2 = (1 - exp(-0.25*(pow(alpha, 2)*log(1/(alpha*GeVfm*r*lambda_bym)) + pow(1 - alpha, 2)*log(1/((1 - alpha)*GeVfm*r*lambda_bym)))*pow(GeVfm*r*Qs, 2)));
  double term3 = (pow(alpha, 2) + pow(1 - alpha, 2))*pow(cyl_bessel_k(1, GeVfm*m*r), 2) + pow(cyl_bessel_k(0, GeVfm*m*r), 2);

  return term1*term2*term3;
}
//__________________________________________________________________________________________


//__________________________________________________________________________________________
//##########################################################################################
//  Find Maximum Probability given quark flavor, random alpha, and Qs
//  Golden Section Search
//##########################################################################################
double Correlator::FindMaximum(double alpha, double m, double Qs, double lower, double upper, double tolerance)
{
  double k = (sqrt(5.) - 1.) / 2.;  //  Golden Ratio
  double xL = upper - k * (upper - lower);  //  Value on the Left
  double xR = lower + k * (upper - lower);  //  Value on the Right

  //  While the difference between the upper and lower bound are greater than the chosen tolerance, keep searching
  while (upper - lower > tolerance)
  {
    // If the value on the left is greater than the value on the right,
    // then make the left position the new Right position
    if (corr(xL, alpha, m, Qs) > corr(xR, alpha, m, Qs))
    {
      upper = xR;
      xR = xL;
      xL = upper - k*(upper - lower);
    }
    // else make the Right position the new Left position
    else
    {
      lower = xL;
      xL = xR;
      xR = lower + k * (upper - lower);
    }
  }
  return corr((lower + upper) / 2., alpha, m, Qs);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Given quark pair data, determine distance and momentum fraction
//##########################################################################################
double Correlator::F(double r, double alpha, double m, double Qs)
{
  //  Since correlation function type was determined in constructor only need to call corr
  return corr(r, alpha, m, Qs);
}
//__________________________________________________________________________________________
