#include "splitting.h"

Splitter::Splitter()
{

}

Splitter::~Splitter()
{

}

Splitter::Splitter(const Splitter &original)
{
  CopySplitter(original);
}

void Splitter::CopySplitter(const Splitter &e)
{
  flavor_chemistry = e.flavor_chemistry;
	dipole_model = e.dipole_model;
	alpha_s = e.alpha_s;
	alpha_min = e.alpha_min;
	r_max = e.r_max;
  e_thresh = e.e_thresh;
  lambda_ = e.lambda_;

}

Splitter& Splitter::operator= (const Splitter& original)
{
	CopySplitter(original);
	return *this;
}

//__________________________________________________________________________________________
//##########################################################################################
//  Select energy of gluon
//##########################################################################################
double Splitter::RollGlue(double e_tot)
{
  double x, y;
  double e_glue = 0;
  bool got_glue = false;
  int num = 0;

  uniform_real_distribution<double> get_energy(e_thresh, e_tot);
  uniform_real_distribution<double> get_probability(0.0, 1/pow(e_thresh, lambda_));

  cout << "begin RollGlue " << get_probability.max() << " " << e_tot << endl;
  while (!got_glue)
  {
    x = get_energy(get_random_number);
    y = get_probability(get_random_number);
    cout << x << " " << y << " " << endl;
    num++;
    if (y < 1/pow(x, lambda_)) got_glue = true;
  }

  return e_glue/e_tot;
}
//__________________________________________________________________________________________


Charge Splitter::RollFlavor(double Qs)
{
  Charge create_charge;

  return create_charge;
}

Quarks Splitter::SplitSample(Sample sampled_energy)
{
  Quarks create_quarks;
  Charge set_charge;
  double gluon_energy_frac;

  cout << "e_tot in SplitSample = " << sampled_energy.e_tot << endl;

  gluon_energy_frac = RollGlue(sampled_energy.e_tot);
  set_charge = RollFlavor(sampled_energy.q_s);

  create_quarks.CreateQuarks(set_charge, gluon_energy_frac, 0.1, 0.1, 0.1);

  return create_quarks;
}
