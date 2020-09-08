#include "correlation.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Correlator
//##########################################################################################
Correlator::Correlator(string model, double lambda)
{
  dipole_model = model;
  lambda_bym = lambda;
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
  dipole_model = e.dipole_model;
  lambda_bym = e.lambda_bym;
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

double Correlator::Vaccum(double r, double alpha, double m, double Qs)
{
  double term1 = r*(pow(m, 2)/(8*pow(M_PI, 2)))*pow(GeVfm, 2);
  double term2 = (1 - exp(-0.25*(pow(alpha, 2) + pow(1 - alpha, 2))*pow(GeVfm*r*Qs, 2)));
  double term3 = (pow(alpha, 2) + pow(1 - alpha, 2))*pow(cyl_bessel_k(1, GeVfm*m*r), 2) + pow(cyl_bessel_k(0, GeVfm*m*r), 2);
  return term1*term2*term3;
}

double Correlator::MVModel(double r, double alpha, double m, double Qs)
{
  double term1 = r*(pow(m, 2)/(8*pow(M_PI, 2)))*pow(GeVfm, 2);
  double term2 = (1 - exp(-0.25*(pow(alpha, 2)*log(1/(alpha*GeVfm*r*lambda_bym)) + pow(1 - alpha, 2)*log(1/((1 - alpha)*GeVfm*r*lambda_bym)))*pow(GeVfm*r*Qs, 2)));
  double term3 = (pow(alpha, 2) + pow(1 - alpha, 2))*pow(cyl_bessel_k(1, GeVfm*m*r), 2) + pow(cyl_bessel_k(0, GeVfm*m*r), 2);

  return term1*term2*term3;
}

function<double(double, double, double, double)> corr;

double Correlator::FindMaximum(double alpha, double m, double Qs, double lower, double upper, double tolerance)
{
  if (dipole_model == "MV")
  {
    corr = bind(&Correlator::MVModel, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4);
  }
  else
  {
    corr = bind(&Correlator::Vaccum, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4);
  }

  double k = (sqrt(5.) - 1.) / 2.;
  double xL = upper - k * (upper - lower);
  double xR = lower + k * (upper - lower);
  while (upper - lower > 0.01)
  {
    if (corr(xL, alpha, m, Qs) > corr(xR, alpha, m, Qs))
    {
      upper = xR;
      xR = xL;
      xL = upper - k*(upper - lower);
    }
    else
    {
      lower = xL;
      xL = xR;
      xR = lower + k * (upper - lower);
    }
  }
  return corr((lower + upper) / 2., alpha, m, Qs);
}

double Correlator::F(double r, double alpha, double m, double Qs)
{
  if (dipole_model == "MV")
  {
    return MVModel(r, alpha, m, Qs);
  }
  else
  {
    return Vaccum(r, alpha, m, Qs);
  }
}
