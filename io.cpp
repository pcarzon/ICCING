#include "io.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Reads in config file and initializes variables
//##########################################################################################
IO::IO(string configFile)
{
  Initialize(); // Initialze config params and map for reading in config params

  //  config input stream
  ifstream input;
  input.open(configFile);

  string var_type;  //  Used to record parameter type read from file

  //******************************************************************************************
  //  Loop through config file
  //******************************************************************************************
  while (!input.eof())
	{
    input >> var_type;  //  Read parameter type

    //  Switch through the possible parameter types
    //  uses var_type as key to map then reads value to class variable
    switch(mapConfigParams[var_type])
    {
      //  case ConfigParam (Does var_type map to ConfigParam?)
      //  input >> input_var (Read in value of var_type)
      //  break; (Stop checking switch and move on)
      case trentoinputdir: input >> trento_input_dir; break;
      case quarkinputfile: input >> quark_input_file; break;
      case eosfile: input >> eos_file; break;
      case outputdir:  input >> output_dir; break;
      case inputtype: input >> input_type; break;
      case outputtype:  input >> output_type; break;

      case eventlabel: input >> event_label; break;
      case firstevent: input >> first_event; break;
      case lastevent: input >> last_event; break;
      case ta:  input >> t_a; break;
      case tb:  input >> t_b; break;

      case kappa:  input >> kappa_; break;
      case rad:  input >> rad_; break;
      case qrad:  input >> qrad_; break;
      case lambda:  input >> lambda_; break;

      case dipolemodel:  input >> dipole_model; break;
      case alphas:  input >> alpha_s; break;
      case alphamin:  input >> alpha_min; break;
      case rmax:  input >> r_max; break;

      case eosemmitlines:  input >> eos_emmit_lines; break;
      case eosscol:  input >> eos_s_col; break;
      case eosecol:  input >> eos_e_col; break;
      case atrento:  input >> a_trento; break;
      case echop:  input >> e_chop; break;

      case bmin:  input >> b_min; break;
      case bmax:  input >> b_max; break;
      case gridmax: input >> grid_max;  break;
      case gridstep:  input >> grid_step; break;
      case tau0:  input >> tau_0; break;
      case ethresh:  input >> e_thresh; break;
      case chargetype:  input >> charge_type; break;

      //#CONFIGPARAM
    }// End of switch
  }// End of while loop

  //  Setting flag for type of charge tracked
  if(charge_type == "BSQ") {  tracked_charge = 0; }
  else if (charge_type == "UDS") {  tracked_charge = 1; }

  current_event = first_event;  //  Set current_event to first_event

  grid_points = 2*(grid_max/grid_step); //  Calculate # grid_points
}// End of Class constructor
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class Destructor
//##########################################################################################
IO::~IO()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  IO Copy Function
//##########################################################################################
void IO::CopyIO(const IO &e)
{
  trento_input_dir = e.trento_input_dir;
  quark_input_file = e.quark_input_file;
  eos_file = e.eos_file;
  output_dir = e.output_dir;
  input_type = e.input_type;
  output_type = e.output_type;

  event_label = e.event_label;
  first_event = e.first_event;
  last_event = e.last_event;
  t_a = e.t_a;
  t_b = e.t_b;

  kappa_ = e.kappa_;
  rad_ = e.rad_;
  qrad_ = e.qrad_;
  lambda_ = e.lambda_;

  dipole_model = e.dipole_model;
  alpha_s = e.alpha_s;
  alpha_min = e.alpha_min;
  r_max = e.r_max;

  eos_emmit_lines = e.eos_emmit_lines;
  eos_s_col = e.eos_s_col;
  eos_e_col = e.eos_e_col;
  a_trento = e.a_trento;
  e_chop = e.e_chop;

  b_min = e.b_min;
  b_max = e.b_max;
  grid_max = e.grid_max;
  grid_step = e.grid_step;
  tau_0 = e.tau_0;
  e_thresh = e.e_thresh;
  charge_type = e.charge_type;
  //#CONFIGPARAM

  tracked_charge = e.tracked_charge;
  current_event = e.current_event;
  grid_points = e.grid_points;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
IO::IO(const IO &original)
{
  CopyIO(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
IO& IO::operator= (const IO& original)
{
	CopyIO(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Initialze variables and map
//##########################################################################################
void IO::Initialize()
{
  //******************************************************************************************
  //  Set variables to default values
  //******************************************************************************************
  trento_input_dir = "";
  quark_input_file = "";
  eos_file = "";
  output_dir = "";
  input_type = 0;
  output_type = 0;

  event_label = "";
  first_event = 0;
  last_event = 0;
  t_a = false;
  t_b = false;

  kappa_ = 0.0;
  rad_ = 0.0;
  qrad_ = 0.0;
  lambda_ = 0.0;

  dipole_model = "";
  alpha_s = 0.0;
  alpha_min = 0.0;
  r_max = 0.0;

  eos_emmit_lines = 0.0;
  eos_s_col = 0;
  eos_e_col = 0;
  a_trento = 0.0;
  e_chop = 0.0;

  b_min = 0;
  b_max = 0;
  grid_max = 0.0;
  grid_step = 0.0;
  tau_0 = 0.0;
  e_thresh = 0.0;
  charge_type = "";
  //#CONFIGPARAM

  tracked_charge = 0;
  current_event = 0;
  grid_points = 0;

  //******************************************************************************************
  //  Initialze map for reading in config file
  //******************************************************************************************
  mapConfigParams["trento_input_dir"] = trentoinputdir;
  mapConfigParams["quark_input_file"] = quarkinputfile;
  mapConfigParams["eos_file"] = eosfile;
  mapConfigParams["output_dir"] = outputdir;
  mapConfigParams["input_type"] = inputtype;
  mapConfigParams["output_type"] = outputtype;

  mapConfigParams["event_label"] = eventlabel;
  mapConfigParams["first_event"] = firstevent;
  mapConfigParams["last_event"] = lastevent;
  mapConfigParams["t_a"] = ta;
  mapConfigParams["t_b"] = tb;

  mapConfigParams["kappa_"] = kappa;
  mapConfigParams["rad_"] = rad;
  mapConfigParams["qrad_"] = qrad;
  mapConfigParams["lambda_"] = lambda;

  mapConfigParams["dipole_model"] = dipolemodel;
  mapConfigParams["alpha_s"] = alphas;
  mapConfigParams["alpha_min"] = alphamin;
  mapConfigParams["r_max"] = rmax;

  mapConfigParams["eos_emmit_lines"] = eosemmitlines;
  mapConfigParams["eos_s_col"] = eosscol;
  mapConfigParams["eos_e_col"] = eosecol;
  mapConfigParams["a_trento"] = atrento;
  mapConfigParams["e_chop"] = echop;

  mapConfigParams["b_min"] = bmin;
  mapConfigParams["b_max"] = bmax;
  mapConfigParams["grid_max"] = gridmax;
  mapConfigParams["grid_step"] = gridstep;
  mapConfigParams["tau_0"] = tau0;
  mapConfigParams["e_thresh"] = ethresh;
  mapConfigParams["charge_type"] = chargetype;
  //#CONFIGPARAM
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Copy config file to ICCING output folder for future reference
//##########################################################################################
void IO::OutputConfig(string file_name)
{
  ofstream output;
  output.open(file_name);

  output
    << "trento_input_dir " << trento_input_dir
    << "\noutput_dir " << output_dir
    << "\ninput_type " << input_type
    << "\noutput_type " << output_type;

  output
    << "\n\nevent_label " << event_label
    << "\nfirst_event " << first_event
    << "\nlast_event " << last_event
    << "\nt_a " << t_a
    << "\nt_b " << t_b;

  output
    << "\n\nkappa_ " << kappa_
    << "\nrad_ " << rad_
    << "\nqrad_ " << qrad_
    << "\nlambda_ " << lambda_;

  output
    << "\n\ndipole_model " << dipole_model
    << "\nalpha_s " << alpha_s
    << "\nalpha_min " << alpha_min
    << "\nr_max " << r_max;

  output
    << "\n\neos_emmit_lines " << eos_emmit_lines
    << "\neos_s_col " << eos_s_col
    << "\neos_e_col " << eos_e_col
    << "\na_trento " << a_trento
    << "\ne_chop " << e_chop;
  output
    << "\n\nb_min " << b_min
    << "\nb_max " << b_max
    << "\ngrid_max " << grid_max
    << "\ngrid_step " << grid_step
    << "\ntau_0 " << tau_0
    << "\ne_thresh " << e_thresh
    << "\ncharge_type " << charge_type;
    //#CONFIGPARAM

    output.close();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Initialize the Event Object
//##########################################################################################
Event IO::InitializeEvent()
{
  Event event_in; //  Temp Event object used to store event specific data
  cout << "We're in!" << endl;
  //  Set variables in event with data from configFile
  event_in.kappa_ = kappa_;
  event_in.lambda_ = lambda_;
  event_in.grid_max = grid_max;
  event_in.grid_step = grid_step;
  event_in.grid_points = grid_points;

  //  Initialize input grid to 0 with dimensions grid_points + 1
  event_in.initial_energy.resize(grid_points + 1, vector<double>(grid_points + 1, 0.));

  event_in.gluon_rad = round(rad_/grid_step);
  event_in.quark_rad = round(qrad_/grid_step);
  cout << event_in.gluon_rad << endl << event_in.quark_rad << endl;

  event_in.gluon_dist.resize(2*event_in.gluon_rad + 3, vector<int>(2*event_in.gluon_rad + 3, 0));
  event_in.quark_dist.resize(2*event_in.quark_rad + 3, vector<double>(2*event_in.quark_rad + 3, 0.));
  cout << event_in.gluon_dist[0].size() << endl;
  for (int i = 0; i < event_in.gluon_dist[0].size(); i++)
  {
    for (int j = event_in.gluon_rad+1; j < event_in.gluon_dist[0].size(); j++)
    {
      //cout << "in again" << endl;
      if (i < event_in.gluon_rad && j < round(sqrt(pow(event_in.gluon_rad,2) - pow(j,2))))
      {
        event_in.gluon_dist[i][j] = 1;
        event_in.gluon_dist[event_in.gluon_rad-i][event_in.gluon_rad-j] = 1;
        cout << event_in.gluon_dist[i][j] << endl;
      }
    }
  }

  ofstream output;
  output.open(output_dir + "gluon_dist_test_1.dat");
  for (int i = 0; i < event_in.gluon_dist[0].size(); i++)
  {
    for (int j = 0; j < event_in.gluon_dist[0].size(); j++)
    {
      output << event_in.gluon_dist[i][j];


      if (j == event_in.gluon_dist[0].size())
      { output << endl;  }  //  If at end of row, go to next row
      else
      { output << " "; }  //  If not at end of row, add space
    }
  }
  output.close();

  return event_in;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Initialize the Splitter Object
//##########################################################################################
Splitter IO::InitializeSplitter()
{
  Splitter init_splitter;

  init_splitter.flavor_chemistry = flavor_chemistry;
  init_splitter.dipole_model = dipole_model;
  init_splitter.alpha_s = alpha_s;
  init_splitter.alpha_min = alpha_min;
  init_splitter.r_max = r_max;

  return init_splitter;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Initialize the ProbDist Object
//##########################################################################################
ProbDist IO::InitializeProbDist()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Initialize the EOS Object
//##########################################################################################
//EOS IO::InitializeEOS()
//{

//}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Print Density Grids with filler 0s
//##########################################################################################
void IO::OutputFullDensityGrids(vector<vector<double>> density_grid, string file_name)
{
  //  Output file stream
  ofstream output;
  output.open(file_name);

  //******************************************************************************************
  //  Loop through density grids
  //******************************************************************************************
  for (int i = 0; i < density_grid.size(); i++)
  {
    for (int j = 0; j < density_grid[0].size(); j++)
    {

      output << density_grid[i][j];  //  Output value at grid point to file

      //  Test to see if at end of row
      if (j == density_grid[0].size() - 1)
      { output << endl;  }  //  If at end of row, go to next row
      else
      { output << " "; }  //  If not at end of row, add space
    }
  }

  output.close();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Print Density Grids without filler 0s
//##########################################################################################
void IO::OutputSparseDensityGrids(vector<vector<double>> density_grid, string file_name)
{
  ofstream output;
  output.open(file_name);

  double x, y, value;

  for (int i = 0; i < density_grid.size(); i++)
  {
    for (int j = 0; j < density_grid[0].size(); j++)
    {
      if (density_grid[i][j] != 0)
      {
        x = -grid_max + i*grid_step;  //  Converts grid point to physical x-value
        y = -grid_max + j*grid_step;  //  Converts grid point to physical y-value
        value = density_grid[i][j];

        output << x << " " << y << " " << value << endl;
      }
    }
  }

  output.close();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Print Eccentricities
//##########################################################################################
void IO::OutputEccentricities(Eccentricity &ecc, string file_name)
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Read Event specific Energy densities
//    Create Event and set densities and variables
//    (Assumes sparse file with only valued points, need to generalize this)
//##########################################################################################
Event IO::ReadEvent()
{
  //  Input file stream
  ifstream input;
  input.open(trento_input_dir + "ic" + to_string(current_event) + ".dat");

  Event event_in; //  Temp Event object used to store event specific data

  // Loop input variables
  int x, y;
  double readx,ready,value;

  //  Ignore first line of input file (Trento specific, needs to be changed)
  input.ignore(10000, '\n');

  //******************************************************************************************
  //  Loop through file until end is reached
  //******************************************************************************************
  while (!input.eof())
  {
      //  Read in point from energy density
      input >> readx >> ready >> value;

      //  Take physical point and convert x and y values into grid indicies
      x = (readx + grid_max)/grid_step;
      y = (ready + grid_max)/grid_step;

      //  Set point in event's initial energy density grid
      event_in.initial_energy[x][y] = value;

      input.ignore(10000, '\n');  //  Ignore rest of line
      if (input.peek() == '\n') {break;}  //  Saftey check for empty line at end of file
  }
  input.close();  //  Close input stream

  //******************************************************************************************
  //  If method requires T_a energy density, read it into event
  //******************************************************************************************
  if (t_a)
  {
    //  Open T_a energy density file
    input.open(trento_input_dir + "TA" + to_string(current_event) + ".dat");
    //  Initialize t_a grid to 0 with dimensions grid_points + 1
    event_in.t_a.resize(grid_points + 1, vector<double>(grid_points + 1, 0));

    //  Ignore first line of input file (Trento specific, needs to be changed)
    input.ignore(10000, '\n');

    //  Loop through file until end is reached
    while (!input.eof())
    {
        //  Read in point from energy density
        input >> readx >> ready >> value;

        //  Take physical point and convert x and y values into grid indicies
        x = (readx + grid_max)/grid_step;
        y = (ready + grid_max)/grid_step;

        //  Set point in event's initial energy density grid
        event_in.t_a[x][y] = value;

        input.ignore(10000, '\n');  //  Ignore rest of line
        if (input.peek() == '\n') {break;}  //  Saftey check for empty line at end of file
    }
    input.close();  //  Close input stream
  }

  //******************************************************************************************
  //  If method requires T_b energy density, read it into event
  //******************************************************************************************
  if (t_b)
  {
    //  Open T_a energy density file
    input.open(trento_input_dir + "TB" + to_string(current_event) + ".dat");
    //  Initialize t_b grid to 0 with dimensions grid_points + 1
    event_in.t_b.resize(grid_points + 1, vector<double>(grid_points + 1, 0));

    //  Ignore first line of input file (Trento specific, needs to be changed)
    input.ignore(10000, '\n');

    //  Loop through file until end is reached
    while (!input.eof())
    {
      //  Read in point from energy density
      input >> readx >> ready >> value;

      //  Take physical point and convert x and y values into grid indicies
      x = (readx + grid_max)/grid_step;
      y = (ready + grid_max)/grid_step;

      //  Set point in event's initial energy density grid
      event_in.t_b[x][y] = value;

      input.ignore(10000, '\n');  //  Ignore rest of line
      if (input.peek() == '\n') {break;}  //  Saftey check for empty line at end of file
    }
    input.close();  //  Close input stream
  }

  return event_in;  //  Return event with data
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Write Energy densities, eccentricities, and configFile
//    Needs to be updated as code is being written
//##########################################################################################
void IO::WriteEvent(Event event)
{
  vector<vector<double>> output_energy;

  //******************************************************************************************
  //  Output Full Density Grids
  //******************************************************************************************
  if (output_type == 0)
  {
    output_energy = event.initial_energy;
    OutputFullDensityGrids(output_energy, output_dir + "ic" + to_string(current_event) + ".dat");

    if (t_a)  //  Output T_a if flag is true
    {
      output_energy = event.t_a;
      OutputFullDensityGrids(output_energy, output_dir + "TA" + to_string(current_event) + ".dat");
    }
    if (t_b)  //  Output T_b if flag is true
    {
      output_energy = event.t_b;
      OutputFullDensityGrids(output_energy, output_dir + "TB" + to_string(current_event) + ".dat");
    }
  }
  //******************************************************************************************
  //  Output Sparse Density Grids
  //******************************************************************************************
  else if (output_type == 1)
  {
    output_energy = event.initial_energy;
    OutputSparseDensityGrids(output_energy, output_dir + "ic" + to_string(current_event) + ".dat");

    if (t_a)  //  Output T_a if flag is true
    {
      output_energy = event.t_a;
      OutputSparseDensityGrids(output_energy, output_dir + "TA" + to_string(current_event) + ".dat");
    }
    if (t_b)  //  Output T_b if flag is true
    {
      output_energy = event.t_b;
      OutputSparseDensityGrids(output_energy, output_dir + "TB" + to_string(current_event) + ".dat");
    }
  }

  current_event++;  //  Used for tracking which event has been processed
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Return true if on last event
//##########################################################################################
bool IO::LastEvent()
{
  if (current_event == last_event)
  { return true;  }
  else
  { return false; }
}
//__________________________________________________________________________________________
