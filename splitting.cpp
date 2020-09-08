#include "splitting.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Splitter
//##########################################################################################
Splitter::Splitter()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
Splitter::~Splitter()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Splitter::Splitter(const Splitter &original)
{
  CopySplitter(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Splitter Copy Function
//##########################################################################################
void Splitter::CopySplitter(const Splitter &e)
{
  flavor_chemistry = e.flavor_chemistry;
	charge_type = e.charge_type;
	alpha_s = e.alpha_s;
	alpha_min = e.alpha_min;
	r_max = e.r_max;
  e_thresh = e.e_thresh;
  lambda_ = e.lambda_;
  grid_step = e.grid_step;
  test_ = e.test_;
  output_dir = e.output_dir;

  Model_Correlator = e.Model_Correlator;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Splitter& Splitter::operator= (const Splitter& original)
{
	CopySplitter(original);
	return *this;
}
//__________________________________________________________________________________________

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

  //  Initialize distributions for selecting gluon energy
  uniform_real_distribution<double> get_energy(e_thresh, e_tot);  //  Selects an energy along the x-axis
  uniform_real_distribution<double> get_probability(0.0, 1.01/pow(e_thresh, lambda_));  //  Selects a probability along the y-axis

  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  Flagged statement for reproducing Gluon Energy Distribution for testing
  if (test_ == "GluonEnergyDist"){ output.open(output_dir + "gluon_energy_dist_test.dat"); }
  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //  Loop until an energy and probability is selected where the probability is less than
  //  the value of the distribution function at selected energy
  while (!got_glue)
  {
    //  Sample random energy and probability
    x = get_energy(get_random_number);
    y = get_probability(get_random_number);

    //  Test if the random probability is less than the function at selected energy
    if (y < 1/pow(x, lambda_)){ got_glue = true; }

    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  Flagged statement for reproducing Gluon Energy Distribution for testing
    if (test_ == "GluonEnergyDist")
    {
      output << x << " " << 1/pow(x, lambda_) << endl;
      num_tests++;
      if (num_tests < 10000){ got_glue = false; }
    }
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  Flagged statement for reproducing Gluon Energy Distribution for testing
  if (test_ == "GluonEnergyDist"){ output.close();  exit(0); }
  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //  Return fraction of energy that is to be selected
  return x/e_tot;
}
//__________________________________________________________________________________________


//__________________________________________________________________________________________
//##########################################################################################
//  Select flavor of gluon
//##########################################################################################
Charge Splitter::RollFlavor(double Qs)
{
  Charge create_charge;
  charge_type = 0;
  //  Initialize distribution for selecting quark flavor probability
  uniform_real_distribution<double> get_flavor(0, 1);

  //  For given Qs, probabilities to get each flavor are extrapolated
  double u = alpha_s*InterpolateValue(FindRange(flavor_chemistry[0], Qs), Qs);
  double d = alpha_s*InterpolateValue(FindRange(flavor_chemistry[1], Qs), Qs);
  double s = alpha_s*InterpolateValue(FindRange(flavor_chemistry[2], Qs), Qs);
  double c = alpha_s*InterpolateValue(FindRange(flavor_chemistry[3], Qs), Qs);

  //  Get probability to get gluon for given Qs
  double g = 1 - u - d - s - c;

  //  Get probability value of gluon to see if it will split
  double probability = get_flavor(get_random_number);

  //  Test if gluon remains gluon
  if (0 <= probability && probability <= g)
  { create_charge.Gluon(charge_type);  }
  //  Test if gluon becomes up quark pair
  else if (g < probability && probability <= g + u)
  { create_charge.Up(charge_type);  }
  //  Test if gluon becomes down quark pair
  else if (g + u < probability && probability <= g + u + d)
  { create_charge.Down(charge_type);  }
  //  Test if gluon becomes strange quark pair
  else if (g + u + d < probability && probability <= g + u + d + s)
  { create_charge.Strange(charge_type);  }
  //  Test if gluon becomes charm quark pair
  else
  { create_charge.Charm(charge_type);  }

  cout << "quark charges " << " " << create_charge.GetCharge()[0] << " " << create_charge.GetCharge()[1] << " " << create_charge.GetCharge()[2] << " " << create_charge.GetCharge()[3] << endl;
  //  Return charge of sample
  return create_charge;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Sample location and momentum fraction of quark pair
//##########################################################################################
vector<double> Splitter::RollLocation(double mass, double Qs)
{
  uniform_real_distribution<double> get_r(0, r_max);
  uniform_real_distribution<double> get_alpha(alpha_min, 1 - alpha_min);
  uniform_real_distribution<double> get_phi(0, 2*M_PI);

  double r, r_final = 0, alpha, phi, prob, ceiling;

  while (r_final == 0)
  {
    r = get_r(get_random_number);
    alpha = get_alpha(get_random_number);

    ceiling = 1.01*Model_Correlator.FindMaximum(alpha, mass, Qs, 0, 1, 0.001);

    uniform_real_distribution<double> get_location_prob(0, ceiling);
    prob = get_location_prob(get_random_number);
    if (prob < Model_Correlator.F(r, alpha, mass, Qs))
    { r_final = r;  }
  }

  phi = get_phi(get_random_number);
  return {alpha, round((r_final*cos(phi))/grid_step), round((r_final*sin(phi))/grid_step)};
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Process energy from Event and create gluon or quark pairs
//##########################################################################################
Quarks Splitter::SplitSample(Sample sampled_energy)
{
  Quarks create_quarks;
  Charge set_charge;
  double gluon_energy_frac;
  vector<double> quark_location;

  //  Get fraction of energy for gluon
  gluon_energy_frac = RollGlue(sampled_energy.e_tot);

  //  Get flavor of gluon
  set_charge = RollFlavor(sampled_energy.q_s);

  //   If this is false go back to SampleEnergy and find new center point
  if (2*set_charge.GetCharge()[0] > gluon_energy_frac*sampled_energy.e_tot)  { gluon_energy_frac = -1; }

  if (set_charge.GetCharge()[0] != 0)
  {
    quark_location = RollLocation(set_charge.GetCharge()[0], sampled_energy.q_s);
  }
  else
  {
    quark_location = {0, 0, 0};
  }

  //  Create quarks to be distributed in output density grids
  create_quarks.CreateQuarks(set_charge, gluon_energy_frac, quark_location[0], quark_location[1], quark_location[2]);

  return create_quarks;
}
//__________________________________________________________________________________________
