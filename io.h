#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>

#include "event.h"
#include "ecc.h"

using namespace std;

class IO
{
private:

  static enum ConfigParams
  {
    inputfile,
    outputfile,
    inputtype,
    outputtype,

    numevents,
    reducedthickness,
    multfluctuations,
    crosssection,
    nucleonwidth,
    bmin,
    bmax,
    gridmax,
    gridstep
  };

  static std::map<string, ConfigParams> mapConfigParams;

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

  static void Initialize();

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
