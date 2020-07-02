#include "io.h"

IO::IO(string configFile)
{

  Initialize();

  ifstream input;
  input.open(configFile);

  string var_type;

  string svar;
  int ivar;
  double dvar;
  string oldvar;

  while (!input.eof())
	{

    input >> oldvar;
    auto var = oldvar;
    switch(mapConfigParams[var_type])
    {
      case inputfile:

        input_file = var;
        break;

      case outputfile:

        output_file = var;
        break;

      case inputtype:

        input_type = var;
        break;

      case outputtype:

        output_type = var;
        break;

      case numevents:

        num_events = var;
        break;

      case reducedthickness:

        reduced_thickness = var;
        break;

      case multfluctuations:

        mult_fluctuations = var;
        break;

      case crosssection:

        cross_section = var;
        break;

      case nucleonwidth:

        nucleon_width = var;
        break;

      case bmin:

        b_min = var;
        break;

      case bmax:

        b_max = var;
        break;

      case gridmax:

        grid_max = var;
        break;

      case gridstep:

        grid_step = var;
        break;
    }
/*    switch(mapConfigParams[var_type])
    {
      case inputfile:
        input >> var;
        input_file = var;
        break;

      case outputfile:
        input >> svar;
        output_file = svar;
        break;

      case inputtype:
        input >> ivar;
        input_type = ivar;
        break;

      case outputtype:
        input >> ivar;
        output_type = ivar;
        break;

      case numevents:
        input >> ivar;
        num_events = ivar;
        break;

      case reducedthickness:
        input >> ivar;
        reduced_thickness = ivar;
        break;

      case multfluctuations:
        input >> dvar;
        mult_fluctuations = dvar;
        break;

      case crosssection:
        input >> dvar;
        cross_section = dvar;
        break;

      case nucleonwidth:
        input >> dvar;
        nucleon_width = dvar;
        break;

      case bmin:
        input >> ivar;
        b_min = ivar;
        break;

      case bmax:
        input >> ivar;
        b_max = ivar;
        break;

      case gridmax:
        input >> dvar;
        grid_max = dvar;
        break;

      case gridstep:
        input >> dvar;
        grid_step = dvar;
        break;
    }*/
  }

  cout << input_file << endl;
  cout << input_type << endl;
  cout << num_events << endl;
  cout << grid_max << endl;
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

Event IO::ReadEvent()
{
  ifstream in;
  in.open(input_file);

  Event inputEvent;

  return inputEvent;
}

void IO::WriteEvent(const Event &event)
{

}
