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
  grid_points = e.grid_points;

  initial_energy = e.initial_energy;
  t_a = e.t_a;
  t_b = e.t_b;
  density = e.density;
  gluon_dist = e.gluon_dist;
  quark_dist = e.quark_dist;

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
  int glue_x_start = 0;
  int glue_y_start = 0;
  int glue_x_end = gluon_dist.size();
  int glue_y_end = gluon_dist.size();
  double q_s = 0;
  double e_tot = 0;

  if (x_center - gluon_rad < 0)
  { glue_x_start = -(x_center - gluon_rad);  }
  if (y_center - gluon_rad < 0)
  { glue_y_start = -(y_center - gluon_rad);  }
//cout << "glue length = " << gluon_dist.size() << endl;
//cout << "glue start x = " << glue_x_start << endl;
//cout << "glue start y = " << glue_y_start << endl;
  if (x_center + gluon_rad > t_b.size())
  { glue_x_end = t_b.size() - (x_center + gluon_rad);  }
  if (y_center + gluon_rad > t_b.size())
  { glue_y_end = t_b.size() - (x_center + gluon_rad);  }
//  cout << "t_b size = " << t_b.size() << endl;
//  cout << "glue end x = " << glue_x_end << endl;
//  cout << "glue end y = " << glue_y_end << endl;

  for (int i = glue_x_start; i < glue_x_end; i++)
  {
    for (int j = glue_y_start; j < glue_y_end; j++)
    {
      q_s += t_b[x_center - gluon_rad + i][y_center - gluon_rad + j]*gluon_dist[i][j];
      e_tot += initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j]*gluon_dist[i][j];
    }
  }
  Sample samp;
  samp.q_s = q_s;
  samp.e_tot = e_tot;
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
  int x = 0;
  int y = 0;
  bool got_point = false;
  int num = 0;
  cout << "begin SampleEnergy " << get_random_number << endl;
  cout << get_grid_point.min() << " " << get_grid_point.max() << endl;
  while (!got_point)
  {
    x = get_grid_point(get_random_number);
    y = get_grid_point(get_random_number);
    cout << x << " " << y << " " << initial_energy[x][y] << endl;
    num++;
    if (initial_energy[x][y] > 0) got_point = true;
  }
  cout << "number of times through loop = " << num << endl;
  cout << x << " " << y << " " << initial_energy[x][y] << endl;

/*  ofstream output;
  output.open("/projects/jnorhos/pcarzon/ICCING/testOutput/gluon_rad_test_0.dat");
  Sample temp;
  for (int i = 0; i < t_b.size()-gluon_rad; i++)
  {
    for (int j = 0; j < t_b.size()-gluon_rad; j++)
    {
      if (initial_energy[i][j] > 0)
      {
        temp = GetQs(i,j);
      }
      if (temp.e_tot > 0)
      {
        output << i*grid_step << " " << j*grid_step << " " << temp.q_s << " " << temp.e_tot << endl;
      }
    }
  }

  output.close();*/
  return GetGlue(x,y);;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Propogates Results of Splitter
//##########################################################################################
void Event::UpdateDensity(Quarks quark_density)
{

}
//__________________________________________________________________________________________

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
