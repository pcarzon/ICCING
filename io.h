#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "event.h"
#include "ecc.h"

using namespace std;

class IO
{
private:

  string input_file;
  string output_file;
  int input_type;
  int output_type;

  int num_events;
  int reduced_thickness;
  double mult_fluctuations;
  double cross_section;
  double nucleon_width;
  int b_min;
  int b_max;
  double grid_max;
  double grid_step;

	void CopyIO(const IO &e);

  void OutputFullDensityGrids(const Event &event);
  void OutputSparseDensityGrids(const Event &event);
  void OutputEccentricities(const Eccentricity &ecc);

public:

  IO(string configFile);
  ~IO();

  IO(const IO &original);
  IO& operator=(const IO& original);

  Event ReadEvent();
  void WriteEvent(const Event &event);
};
#endif
