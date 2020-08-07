#include "event.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Event
//##########################################################################################
Event::Event()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
Event::~Event()
{
  CleanEvent();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Event::Event(const Event &original)
{
	CopyEvent(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Event Copy Function
//##########################################################################################
void Event::CopyEvent(const Event &e)
{
//  cout << "started copying event" << endl;

  kappa_ = e.kappa_;
  gluon_rad = e.gluon_rad;
  quark_rad = e.quark_rad;
  lambda_ = e.lambda_;
  grid_max = e.grid_max;
  grid_step = e.grid_step;
  tau_0 = e.tau_0;
  e_thresh = e.e_thresh;
  grid_points = e.grid_points;
  get_grid_point = e.get_grid_point;

  initial_energy = e.initial_energy;
  t_a = e.t_a;
  t_b = e.t_b;
  density = e.density;
  gluon_dist = e.gluon_dist;
  quark_dist = e.quark_dist;
  valued_points = e.valued_points;

  total_initial_energy = e.total_initial_energy;
  total_energy = e.total_energy;
  seed = e.seed;

//  cout << "Finished copying event" << endl;

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Event& Event::operator= (const Event& original)
{
	CopyEvent(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Select energy of gluon
//##########################################################################################
Sample Event::GetGlue(int x_center, int y_center)
{
  int total_points = 0;
  vector<int> gluon_bounds = GetIntegrationBounds(x_center, y_center, gluon_rad);

  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      q_s += t_b[x_center - gluon_rad + i][y_center - gluon_rad + j]*gluon_dist[i][j];
      e_tot += initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j]*gluon_dist[i][j];
      total_points++;
    }
  }

  Sample samp;
  samp.q_s = q_s/total_points;
  samp.e_tot = pow(grid_step,2)*tau_0*e_tot;

  return samp;
}
//__________________________________________________________________________________________


//__________________________________________________________________________________________
//##########################################################################################
//  Sample Initial Energy for ICCING algorithm
//    See: First 2 commands in While in DistributeCharge in ICCING_v0_1_8.nb, Calls RollGlue
//##########################################################################################
Sample Event::SampleEnergy()
{
  int point;
  bool got_point = false;
  int num = 0;
  Sample out_sample;

  while (!got_point)
  {
    //  Initialize distribution for selecting points from initial_energy
    get_grid_point = uniform_int_distribution<int>(0, valued_points.size());

    point = get_grid_point(get_random_number);
    num++;

    out_sample = GetGlue(valued_points[point][0], valued_points[point][1]);

  //    if (out_sample.e_tot < e_thresh)
  //    {
        UpdateEnergy(valued_points[point][0], valued_points[point][1], 1.);
  //      continue;
  //    }
  //    else
  //    {
        got_point = true;
  //    }
  }
//  cout << "number of times through loop = " << num << endl;
//  cout << valued_points[point][0] << " " << valued_points[point][1] << " " << initial_energy[valued_points[point][0]][valued_points[point][1]] << endl;

  return out_sample;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Propogates Results of Splitter
//##########################################################################################
void Event::UpdateDensity(Quarks quark_density)
{
  vector<double> position = quark_density.GetPosition()
  if (quark_density.GetCharge()[0] == 0.)
  {
    UpdateEnergy(position[0], position[1], quark_density.GetEnergyFraction());
  }

  vector<int> gluon_bounds = GetIntegrationBounds(x_center, y_center, gluon_rad);

  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      initial_energy[i][j] -= gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
    }
  }

  vector<int> quark_bounds = GetIntegrationBounds(x_center, y_center, quark_rad);

  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      density[1][i][j] -= quark_dist[i][j]*initial_energy[x_center - quark_rad + i][y_center - quark_rad + j];
      density[2][i][j] -= quark_dist[i][j]*initial_energy[x_center - quark_rad + i][y_center - quark_rad + j];

      if (quark_density.GetCharge()[2] == -1)
      {
        density[3][i][j] -= quark_dist[i][j]*initial_energy[x_center - quark_rad + i][y_center - quark_rad + j];
      }
  }
  }

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Subtracts energy from initial_energy and adds it to density[0]
//##########################################################################################
void Event::UpdateEnergy(int x_center, int y_center, double ratio)
{
  vector<int> gluon_bounds = GetIntegrationBounds(x_center, y_center, gluon_rad);

  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      initial_energy[i][j] -= gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
      density[0][i][j] += gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
    }
  }
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Gets intigration bounds for density grid manipulations
//##########################################################################################
vector<int> Event::GetIntegrationBounds(int x_center, int y_center, double raduis);
{
  int glue_x_start = 0;
  int glue_y_start = 0;
  int glue_x_end = gluon_dist.size();
  int glue_y_end = gluon_dist.size();
  vector<int> bounds;
  bounds.push_back(0);
  bounds.push_back(grid_points);
  bounds.push_back(0);
  bounds.push_back(grid_points);

  if (x_center - gluon_rad < 0)
  { bounds[0] = -(x_center - gluon_rad);  }
  if (y_center - gluon_rad < 0)
  { bounds[2] = -(y_center - gluon_rad);  }

  if (x_center + gluon_rad > t_b.size())
  { bounds[1] = t_b.size() - (x_center + gluon_rad);  }
  if (y_center + gluon_rad > t_b.size())
  { bounds[3] = t_b.size() - (x_center + gluon_rad);  }

  return bounds;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Checks Event totals and returns true when initial_total is below a threshold
//##########################################################################################
bool Event::IsEventDone()
{
  for (int i = 0; i < valued_points.size(); i++)
  {
    if (initial_energy[valued_points[i][0]][valued_points[i][1]] == 0)
    { valued_points.erase(valued_points.begin() + i); }
  }
  if (valued_points.size() == 0)
  { return true;  }
  else
  { return false; }
}

//__________________________________________________________________________________________
//##########################################################################################
//  Clears Event variables as a cautionary measure
//##########################################################################################
void Event::CleanEvent()
{
  initial_energy.clear();
  t_a.clear();
  t_b.clear();
  density.clear();
  gluon_dist.clear();
  quark_dist.clear();

  total_initial_energy = 0;
  total_energy = 0;
  seed = 0;
}
//__________________________________________________________________________________________
