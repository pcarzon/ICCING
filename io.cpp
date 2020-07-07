#include "io.h"

//  Class constructor, must specify path to configuration file
IO::IO(string configFile)
{
  // Initialze map for reading in config parameters
  Initialize();
cout << "Passed Initialize\n";
  //  config input stream
  ifstream input;
  input.open(configFile);
  cout << "Opened File\n";

  //  Used to record parameter type read from file
  string var_type;

  //  Loop through config file
  while (!input.eof())
	{
    //  Read parameter type
    input >> var_type;
    cout << "Reading Config\n";

    //  Switch through the possible parameter types
    //  uses var_type as key to map then reads value to class variable
    switch(mapConfigParams[var_type])
    {
      case inputfolder: input >> input_folder; break;

      case outputfolder:  input >> output_folder; break;

      case inputtype: input >> input_type; break;

      case outputtype:  input >> output_type; break;

      case firstevent: input >> first_event; break;

      case numevents: input >> num_events; break;

      case reducedthickness:  input >> reduced_thickness; break;

      case multfluctuations:  input >> mult_fluctuations; break;

      case crosssection:  input >> cross_section; break;

      case nucleonwidth:  input >> nucleon_width; break;

      case bmin:  input >> b_min; break;

      case bmax:  input >> b_max; break;

      case gridmax: input >> grid_max;  break;

      case gridstep:  input >> grid_step; break;

      case ta:  input >> t_a; cout << "Read t_a " << t_a << endl;
 break;

      case tb:  input >> t_b; break;
      //#CONFIGPARAM
    }// End of switch
  }// End of while loop
  cout << "Read Config\n";

  current_event = first_event;
  //  Calculate # grid_points from config file
  grid_points = 2*(grid_max/grid_step);
}

IO::~IO()
{
  input_folder = "";
  output_folder = "";
  input_type = 0;
  output_type = 0;

  first_event = 0;
  num_events = 0;
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
}

IO::IO(const IO &original)
{
  CopyIO(original);
}

void IO::CopyIO(const IO &e)
{
  input_folder = e.input_folder;
  output_folder = e.output_folder;
  input_type = e.input_type;
  output_type = e.output_type;

  first_event = e.first_event;
  num_events = e.num_events;
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

IO& IO::operator= (const IO& original)
{
	CopyIO(original);
	return *this;
}

void IO::Initialize()
{
  input_folder = "";
  output_folder = "";
  input_type = 0;
  output_type = 0;

  first_event = 0;
  num_events = 0;
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

  mapConfigParams["input_folder"] = inputfolder;
  mapConfigParams["output_folder"] = outputfolder;
  mapConfigParams["input_type"] = inputtype;
  mapConfigParams["output_type"] = outputtype;
  mapConfigParams["first_event"] = firstevent;
  mapConfigParams["num_events"] = numevents;
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

void IO::OutputConfig(string file_name)
{
  ofstream output;
  output.open(file_name);

  output << "input_folder " << input_folder
    << "\noutput_folder " << output_folder
    << "\ninput_type " << input_type
    << "\noutput_type " << output_type;

  output << "\n\nfirst_event " << first_event
    << "\n\nnum_events " << num_events
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

void IO::OutputFullDensityGrids(vector<vector<double>> density_grid, string file_name)
{
  ofstream output;
  output.open(file_name);

  double x, y, value;

  for (int i = 0; i < density_grid.size(); i++)
  {
    for (int j = 0; j < density_grid[0].size(); j++)
    {
      x = -grid_max + i*grid_step;
      y = -grid_max + j*grid_step;
      value = density_grid[i][j];

      output << x << " " << y << " " << value;

      if (j == density_grid[0].size() - 1)
      { output << endl;  }
      else
      { output << " "; }
    }
  }

  output.close();
}

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
        x = -grid_max + i*grid_step;
        y = -grid_max + j*grid_step;
        value = density_grid[i][j];

        output << x << " " << y << " " << value << endl;
      //  cout << x << " " << y << " " << value << endl;
      }
    }
  }

  output.close();
}

void IO::OutputEccentricities(Eccentricity &ecc, string file_name)
{

}

Event IO::ReadEvent()
{
  ifstream input;
  input.open(input_folder + "ic" + to_string(current_event) + ".dat");
  cout << input_folder + "ic" + to_string(current_event) + ".dat" << endl;
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
    input.open(input_folder + "TA" + to_string(current_event) + ".dat");
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
    input.open(input_folder + "TB" + to_string(current_event) + ".dat");
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

  current_event++;

  return event_in;
}

void IO::WriteEvent(Event &event)
{
  vector<vector<double>> output_energy;

  if (output_type == 0)
  {
    output_energy = event.GetInitialEnergy();
    OutputFullDensityGrids(output_energy, output_folder + "ic" + to_string(current_event) + ".dat");

    if (t_a)
    {
      output_energy = event.GetTa();
      OutputFullDensityGrids(output_energy, output_folder + "TA" + to_string(current_event) + ".dat");
    }
    if (t_b)
    {
      output_energy = event.GetTb();
      OutputFullDensityGrids(output_energy, output_folder + "TB" + to_string(current_event) + ".dat");
    }
  }
  else if (output_type == 1)
  {
    output_energy = event.GetInitialEnergy();
    OutputSparseDensityGrids(output_energy, output_folder + "ic" + to_string(current_event) + ".dat");

    if (t_a)
    {
      output_energy = event.GetTa();
      OutputSparseDensityGrids(output_energy, output_folder + "TA" + to_string(current_event) + ".dat");
    }
    if (t_b)
    {
      output_energy = event.GetTb();
      OutputSparseDensityGrids(output_energy, output_folder + "TB" + to_string(current_event) + ".dat");
    }
  }

}

  bool IO::LastEvent()
  {
    if (current_event == num_events)
    { return true;  }
    else
    { return false; }
  }
