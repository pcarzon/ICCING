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
  test_ = e.test_;
  output_dir = e.output_dir;
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
  bool got_glue = false;
  int num_tests = 0;
  ofstream output;

  uniform_real_distribution<double> get_energy(e_thresh, e_tot);
  uniform_real_distribution<double> get_probability(0.0, 1.01/pow(e_thresh, lambda_));
  cout << "got here " << test_ << endl;
  if (test_ == "GluonEnergyDist")
  {
    cout << "test works " << test_ << endl;
    output.open(output_dir + "gluon_energy_dist_test.dat");

    output << "max prob: " << get_probability.max() << " e_tot: " << e_tot << endl;
  }

  while (!got_glue)
  {
    x = get_energy(get_random_number);
    y = get_probability(get_random_number);

    if (y < 1/pow(x, lambda_)) got_glue = true;// output << "\tyes" << endl;}

    if (test_ == "GluonEnergyDist")
    {
      output << x << " " << 1/pow(x, lambda_) << endl;

      num_tests++;
      if (num_tests < 10000)
      { got_glue = false; }
    }
}

  if (test_ == "GluonEnergyDist")
  {
    output.close();
    exit(0);
  }

  return x/e_tot;
}
//__________________________________________________________________________________________


Charge Splitter::RollFlavor(double Qs)
{
  Charge create_charge;

  uniform_real_distribution<double> get_flavor(0, 1);

  double u = alpha_s*InterpolateValue(FindRange(flavor_chemistry[0], Qs), Qs);
  double d = alpha_s*InterpolateValue(FindRange(flavor_chemistry[1], Qs), Qs);
  double s = alpha_s*InterpolateValue(FindRange(flavor_chemistry[2], Qs), Qs);
  double c = alpha_s*InterpolateValue(FindRange(flavor_chemistry[3], Qs), Qs);
  double g = 1 - u - d - s - c;
  double probability = get_flavor(get_random_number);
  probability = 0.985;
  // gluon prob = 1 - sum of q_s_range
  if (0 <= probability && probability <= g)
  { create_charge.Gluon(dipole_model);  }
  else if (g < probability && probability <= g + u)
  { create_charge.Up(dipole_model);  }
  else if (g + u < probability && probability <= g + u + d)
  { create_charge.Down(dipole_model);  }
  else if (g + u + d < probability && probability <= g + u + d + s)
  { create_charge.Strange(dipole_model);  }
  else
  { create_charge.Charm(dipole_model);  }

  cout << "qs: "
  << FindRange(flavor_chemistry[0], Qs).x << " "
  << FindRange(flavor_chemistry[1], Qs).x << " "
  << FindRange(flavor_chemistry[2], Qs).x << " "
  << FindRange(flavor_chemistry[3], Qs).x << endl;

  cout << "flavor probs: "
  << probability << " "
  << Qs << " "
  << g << " "
  << u << " "
  << d << " "
  << s << " "
  << c << endl;
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
  // If statement to check if 2*quark_mass < gluon_energy_frac*e_tot (add now)
  //    If this is false go back to SampleEnergy and find new center point
  cout << "Flavor " << set_charge.GetCharge()[0] << endl;
  create_quarks.CreateQuarks(set_charge, gluon_energy_frac, 0.1, 0.1, 0.1);

  return create_quarks;
}
