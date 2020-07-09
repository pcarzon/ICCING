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
      case outputdir:  input >> output_dir; break;
      case inputtype: input >> input_type; break;
      case outputtype:  input >> output_type; break;
      case firstevent: input >> first_event; break;
      case lastevent: input >> last_event; break;
      case reducedthickness:  input >> reduced_thickness; break;
      case multfluctuations:  input >> mult_fluctuations; break;
      case crosssection:  input >> cross_section; break;
      case nucleonwidth:  input >> nucleon_width; break;
      case bmin:  input >> b_min; break;
      case bmax:  input >> b_max; break;
      case gridmax: input >> grid_max;  break;
      case gridstep:  input >> grid_step; break;
      case ta:  input >> t_a; break;
      case tb:  input >> t_b; break;
      //#CONFIGPARAM
    }// End of switch
  }// End of while loop

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
  output_dir = e.output_dir;
  input_type = e.input_type;
  output_type = e.output_type;

  first_event = e.first_event;
  last_event = e.last_event;
  reduced_thickness = e.reduced_thickness;
  mult_fluctuations = e.mult_fluctuations;
  cross_section = e.cross_section;
  nucleon_width = e.nucleon_width;
  b_min = e.b_min;
  b_max = e.b_max;
  grid_max = e.grid_max;
  grid_step = e.grid_step;

  t_a = e.t_a;
  t_b = e.t_b;
  //#CONFIGPARAM
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
  output_dir = "";
  input_type = 0;
  output_type = 0;

  first_event = 0;
  last_event = 0;
  reduced_thickness = 0;
  mult_fluctuations = 0.0;
  cross_section = 0.0;
  nucleon_width = 0.0;
  b_min = 0;
  b_max = 0;
  grid_max = 0.0;
  grid_step = 0.0;

  t_a = false;
  t_b = false;
  //#CONFIGPARAM

  //******************************************************************************************
  //  Initialze map for reading in config file
  //******************************************************************************************
  mapConfigParams["trento_input_dir"] = trentoinputdir;
  mapConfigParams["output_dir"] = outputdir;
  mapConfigParams["input_type"] = inputtype;
  mapConfigParams["output_type"] = outputtype;
  mapConfigParams["first_event"] = firstevent;
  mapConfigParams["last_event"] = lastevent;
  mapConfigParams["reduced_thickness"] = reducedthickness;
  mapConfigParams["mult_fluctuations"] = multfluctuations;
  mapConfigParams["cross_section"] = crosssection;
  mapConfigParams["nucleon_width"] = nucleonwidth;
  mapConfigParams["b_min"] = bmin;
  mapConfigParams["b_max"] = bmax;
  mapConfigParams["grid_max"] = gridmax;
  mapConfigParams["grid_step"] = gridstep;
  mapConfigParams["t_a"] = ta;
  mapConfigParams["t_b"] = tb;
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

  output << "trento_input_dir " << trento_input_dir
    << "\noutput_dir " << output_dir
    << "\ninput_type " << input_type
    << "\noutput_type " << output_type;

  output << "\n\nfirst_event " << first_event
    << "\n\nlast_event " << last_event
    << "\nreduced_thickness " << reduced_thickness
    << "\nmult_fluctuations " << mult_fluctuations
    << "\ncross_section " << cross_section
    << "\nnucleon_width " << nucleon_width
    << "\nb_min " << b_min
    << "\nb_max " << b_max
    << "\ngrid_max " << grid_max
    << "\ngrid_step " << grid_step;

    output << "\n\nt_a " << t_a
      << "\nt_b " << t_b;
    //#CONFIGPARAM

    output.close();
}
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

void IO::OutputEccentricities(Eccentricity &ecc, string file_name)
{

}

Event IO::ReadEvent()
{
  ifstream input;
  input.open(trento_input_dir + "ic" + to_string(current_event) + ".dat");
  cout << trento_input_dir + "ic" + to_string(current_event) + ".dat" << endl;
  Event event_in;

  event_in.grid_max = grid_max;
  event_in.grid_step = grid_step;
  event_in.grid_points = grid_points;

  //vector<vector<double>> carrier;
  //  Initialize input grid to 0 with dimensions grid_points
  event_in.initial_energy.resize(grid_points + 1, vector<double>(grid_points + 1, 0));

  int x, y;
  double readx,ready,value;

  input.ignore(10000, '\n');

  while (!input.eof())
  {

      input >> readx >> ready >> value;

      x = (readx + grid_max)/grid_step;
      y = (ready + grid_max)/grid_step;
      //cout << x << " " << y << " " << value << endl;
      event_in.initial_energy[x][y] = value;

      input.ignore(10000, '\n');
      if (input.peek() == '\n') {break;}
  }
  input.close();

  if (t_a)
  {
    input.open(trento_input_dir + "TA" + to_string(current_event) + ".dat");
    event_in.t_a.resize(grid_points + 1, vector<double>(grid_points + 1, 0));

    input.ignore(10000, '\n');

    while (!input.eof())
    {

        input >> readx >> ready >> value;

        x = (readx + grid_max)/grid_step;
        y = (ready + grid_max)/grid_step;
        //cout << x << " " << y << " " << value << endl;
        event_in.t_a[x][y] = value;

        input.ignore(10000, '\n');
        if (input.peek() == '\n') {break;}
    }
    input.close();
  }
  if (t_b)
  {
    input.open(trento_input_dir + "TB" + to_string(current_event) + ".dat");
    event_in.t_b.resize(grid_points + 1, vector<double>(grid_points + 1, 0));

    input.ignore(10000, '\n');

    while (!input.eof())
    {

        input >> readx >> ready >> value;

        x = (readx + grid_max)/grid_step;
        y = (ready + grid_max)/grid_step;
        //cout << x << " " << y << " " << value << endl;
        event_in.t_b[x][y] = value;

        input.ignore(10000, '\n');
        if (input.peek() == '\n') {break;}
    }
    input.close();
  }

  return event_in;
}

void IO::WriteEvent(Event event)
{
  vector<vector<double>> output_energy;

  if (output_type == 0)
  {
    output_energy = event.GetInitialEnergy();
    OutputFullDensityGrids(output_energy, output_dir + "ic" + to_string(current_event) + ".dat");

    if (t_a)
    {
      output_energy = event.GetTa();
      OutputFullDensityGrids(output_energy, output_dir + "TA" + to_string(current_event) + ".dat");
    }
    if (t_b)
    {
      output_energy = event.GetTb();
      OutputFullDensityGrids(output_energy, output_dir + "TB" + to_string(current_event) + ".dat");
    }
  }
  else if (output_type == 1)
  {
    output_energy = event.GetInitialEnergy();
    OutputSparseDensityGrids(output_energy, output_dir + "ic" + to_string(current_event) + ".dat");

    if (t_a)
    {
      output_energy = event.GetTa();
      OutputSparseDensityGrids(output_energy, output_dir + "TA" + to_string(current_event) + ".dat");
    }
    if (t_b)
    {
      output_energy = event.GetTb();
      OutputSparseDensityGrids(output_energy, output_dir + "TB" + to_string(current_event) + ".dat");
    }
  }

  current_event++;
}

bool IO::LastEvent()
{
  if (current_event == last_event)
  { return true;  }
  else
  { return false; }
}
