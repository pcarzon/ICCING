#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include "event.h"
using namespace std;

//class Event;
class Eccentricity;

class IO
{
private:

  //  This is equivalent to a data type with values that can be those inside
  //  ConfigParams is used to identify parameters from the config file
  //    TO ADD CONFIG PARAMETER SEARCH THIS TAG: #CONFIGPARAM
  enum ConfigParams
  {
    inputfolder,
    outputfolder,
    inputtype,
    outputtype,

    firstevent,
    numevents,
    reducedthickness,
    multfluctuations,
    crosssection,
    nucleonwidth,
    bmin,
    bmax,
    gridmax,
    gridstep,

    ta,
    tb
    //#CONFIGPARAM
  };

  //  This maps the values of ConfigParams to keyword strings used in the config
  //  file.
  std::map<string, ConfigParams> mapConfigParams;

  //  These are the file variables
  string input_folder;
  string output_folder;
  int input_type;
  int output_type;

  //  These are variables that are specified in the configuration file
  int first_event;
  int current_event;
  int num_events;
  int reduced_thickness;
  double mult_fluctuations;
  double cross_section;
  double nucleon_width;
  int b_min;
  int b_max;
  double grid_max;
  double grid_step;

  bool t_a;
  bool t_b;
  //#CONFIGPARAM

  // This is calculated to be 2*(grid_max/grid_step)
  int grid_points;

  Event* event_in;

  //  This function initializes the map used to read the config file
  void Initialize();

  //  Copy function for IO class, called by the operator= and implicit copy functions
	void CopyIO(const IO &e);

  //  These functions determine the possible output formats
  void OutputConfig(string file_name);
  void OutputFullDensityGrids(vector<vector<double>> density_grid, string file_name);
  void OutputSparseDensityGrids(vector<vector<double>> density_grid, string file_name);  //  Only prints valued points
  void OutputEccentricities(Eccentricity &ecc, string file_name);

public:

  //  Class constructor, must specify path to configuration file
  IO(string configFile);
  //  Destructor, used to clear all data stored by class
  ~IO();

  //  Implicit copy function, newIOObject(oldIOObject)
  IO(const IO &original);
  //  This function defines what happens when you use the = operator on this class
  IO& operator=(const IO& original);

  //  Reads a single event
  Event* ReadEvent();
  //  Writes a single event to output file
  void WriteEvent(Event* event);

  bool LastEvent();
};
#endif
