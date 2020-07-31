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

  uniform_real_distribution<double> get_flavor(0, 1);

  vector<double> q_s_range;
  for (int i = 0; i < 4; i++)
  {
    q_s_range.push_back(alpha_s*InterpolateValue(FindRange(flavor_chemistry[i], Qs), Qs));
  }

  double probability = get_flavor(get_random_number);
  vector<double>::iterator q_s_prob = q_s_range.begin();

  if (0 <= probability && probability <= *q_s_prob)
  { create_charge.Gluon(dipole_model);  }
  else if (*q_s_prob < probability && probability <= accumulate(q_s_prob, q_s_prob + 1, 0))
  { create_charge.Up(dipole_model);  }
  else if (accumulate(q_s_prob, q_s_prob + 1, 0) < probability && probability <= accumulate(q_s_prob, q_s_prob + 2, 0))
  { create_charge.Down(dipole_model);  }
  else if (accumulate(q_s_prob, q_s_prob + 2, 0) < probability && probability <= accumulate(q_s_prob, q_s_prob + 3, 0))
  { create_charge.Strange(dipole_model);  }
  else
  { create_charge.Charm(dipole_model);  }
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

  cout << "Flavor " << set_charge.GetCharge()[0] << endl;
  create_quarks.CreateQuarks(set_charge, gluon_energy_frac, 0.1, 0.1, 0.1);

  return create_quarks;
}
