#include "io.h"

IO::IO(string configFile)
{

  ifstream input;
  input.open(teamFile);

  string var_type;

  string svar;
  int ivar;
  double dvar;

  while (!input.eof())
	{

    input >> var_type;

    switch(var_type)
    {
      case "input_file":
        input >> svar;
        input_file = svar;
        break;

      case "output_file":
        input >> svar;
        output_file = svar;
        break;

      case "input_type":
        input >> ivar;
        input_type = ivar;
        break;

      case "output_type":
        input >> ivar;
        output_type = ivar;
        break;

      case "num_events":
        input >> ivar;
        num_events = ivar;
        break;

      case "reduced_thickness":
        input >> ivar;
        reduced_thickness = ivar;
        break;

      case "mult_fluctuations":
        input >> dvar;
        mult_fluctuations = dvar;
        break;

      case "cross_section":
        input >> dvar;
        cross_section = dvar;
        break;

      case "nucleon_width":
        input >> dvar;
        nucleon_width = dvar;
        break;

      case "b_min":
        input >> ivar;
        b_min = ivar;
        break;

      case "b_max":
        input >> ivar;
        b_max = ivar;
        break;

      case "grid_max":
        input >> dvar;
        grid_max = dvar;
        break;

      case "grid_step":
        input >> dvar;
        grid_step = dvar;
        break;
    }
  }
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
