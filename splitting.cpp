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
  double e_glue = e_tot;

  return e_glue;
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

  set_charge = RollFlavor(sampled_energy.q_s);

  create_quarks.CreateQuarks(set_charge, 0.1, 0.1, 0.1);

  return create_quarks;
}
