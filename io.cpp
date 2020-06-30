#include "io.h"

IO::IO(string input, string output, int iType, int oType)
{
  input_file = input;
  output_file = output;
  input_type = iType;
  output_type = oType;
}

IO::~IO()
{
  input_file = "";
  output_file = "";
  input_type = 0;
  output_type = 0;
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
