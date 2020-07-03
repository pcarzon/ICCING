#include "io.h"

IO::IO(string configFile)
{

  Initialize();

  ifstream input;
  input.open(configFile);

  string var_type;

  while (!input.eof())
	{
    input >> var_type;
    switch(mapConfigParams[var_type])
    {
      case inputfile:
        input >> input_file; break;

      case outputfile:
        input >> output_file; break;

      case inputtype:
        input >> input_type; break;

      case outputtype:
        input >> output_type; break;

      case numevents:
        input >> num_events; break;

      case reducedthickness:
        input >> reduced_thickness; break;

      case multfluctuations:
        input >> mult_fluctuations; break;

      case crosssection:
        input >> cross_section; break;

      case nucleonwidth:
        input >> nucleon_width; break;

      case bmin:
        input >> b_min; break;

      case bmax:
        input >> b_max; break;

      case gridmax:
        input >> grid_max;  break;

      case gridstep:
        input >> grid_step; break;
    }
  }

  grid_points = 2*(grid_max/grid_step);
  carrier.resize(grid_points, vector<double>(grid_points, 0));
  cout << "carrier size " << carrier.size() << " " << carrier[0].size() << endl;
}

IO::~IO()
{
  input_file = "";
  output_file = "";
  input_type = 0;
  output_type = 0;

  num_events = 0;
  reduced_thickness = 0;
  mult_fluctuations = 0.0;
  cross_section = 0.0;
  nucleon_width = 0.0;
  b_min = 0;
  b_max = 0;
  grid_max = 0.0;
  grid_step = 0.0;
}

IO::IO(const IO &original)
{
  CopyIO(original);
}

void IO::CopyIO(const IO &e)
{
  input_file = e.input_file;
  output_file = e.output_file;
  input_type = e.input_type;
  output_type = e.output_type;

  num_events = e.num_events;
  reduced_thickness = e.reduced_thickness;
  mult_fluctuations = e.mult_fluctuations;
  cross_section = e.cross_section;
  nucleon_width = e.nucleon_width;
  b_min = e.b_min;
  b_max = e.b_max;
  grid_max = e.grid_max;
  grid_step = e.grid_step;
}

IO& IO::operator= (const IO& original)
{
	CopyIO(original);
	return *this;
}

void IO::Initialize()
{
  mapConfigParams["input_file"] = inputfile;
  mapConfigParams["output_file"] = outputfile;
  mapConfigParams["input_type"] = inputtype;
  mapConfigParams["output_type"] = outputtype;
  mapConfigParams["num_events"] = numevents;
  mapConfigParams["reduced_thickness"] = reducedthickness;
  mapConfigParams["mult_fluctuations"] = multfluctuations;
  mapConfigParams["cross_section"] = crosssection;
  mapConfigParams["nucleon_width"] = nucleonwidth;
  mapConfigParams["b_min"] = bmin;
  mapConfigParams["b_max"] = bmax;
  mapConfigParams["grid_max"] = gridmax;
  mapConfigParams["grid_step"] = gridstep;
}

void IO::OutputFullDensityGrids(const Event &event)
{

}

void IO::OutputSparseDensityGrids(const Event &event)
{

}

void IO::OutputEccentricities(const Eccentricity &ecc)
{

}

vector<vector<double>> IO::ReadEvent()
{
  ifstream input;
  input.open(input_file);

  int x, y;
  double value;
	input.ignore(10000, '\n');

  while (!input.eof())
  {
      cout << "Here it is " << int(input.peek()) << endl;
      input >> value;
      x = value/grid_step - 1 + grid_points/2;
      cout << x << endl;
      input >> value;

      y = value/grid_step - 1 + grid_points/2;
      cout << y << endl;
      input >> value;
      carrier[x][y] = value;
      cout << value << carrier[x][y] << endl;
  }
  return carrier;
}

void IO::WriteEvent(const Event &event)
{

}
