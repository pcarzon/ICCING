#ifndef IO_H
#define IO_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <random>
#include <ctime>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################
#include "event.h"
#include "ecc.h"
#include "splitting.h"
#include "functions.h"
#include "correlation.h"
//__________________________________________________________________________________________

using namespace std;
extern default_random_engine get_random_number;

class IO
{
private:
//__________________________________________________________________________________________
//##########################################################################################
//  Configuration Parameters
//##########################################################################################
//    TO ADD CONFIG PARAMETER SEARCH THIS TAG: #CONFIGPARAM
  enum ConfigParams
  {
    //  Equivalent to a data type with values that can be those inside
    //  ConfigParams is used to identify parameters from config file
    trentoinputdir,
    quarkinputfile,
    eosfile,
    outputdir,
    inputtype,
    outputtype,
    seed,
    test,

    eventlabel,
    firstevent,
    lastevent,
    ta,
    tb,

    kappa,
    rad,
    qrad,
    lambda,

    dipolemodel,
    alphas,
    alphamin,
    rmax,
    lambdabym,

    atrento,
    schop,

    gridmax,
    gridstep,
    tau0,
    ethresh,
    chargetype
    //#CONFIGPARAM
  };

  //  Maps values of ConfigParams to keyword strings used in config file
  map<string, ConfigParams> mapConfigParams;

  //  File variables
  string trento_input_dir;
  string quark_input_file;
  string eos_file;
  string output_dir;
  int input_type; //  Type of input: 0 = Full Density Grid, 1 = Sparse Density Grid
  int output_type;  // Type of output: 0 = Full Density Grids, 1 = Sparse Density Grids
  int seed_;
  string test_;

  //******************************************************************************************
  //  Internal Config parameters
  //******************************************************************************************
  string event_label;
  int first_event;  //  First event to be read i.e. 0
  int last_event; //  Last event to be read
  bool t_a; // Is T_a read in
  bool t_b; // Is T_b read in

  //******************************************************************************************
  //  Event Config parameters
  //******************************************************************************************
  double kappa_; //  Used for Qs grid
  double rad_;
  double qrad_;
  double lambda_;

  //******************************************************************************************
  //  Splitter Config parameters
  //******************************************************************************************
  string dipole_model;  //  Future Update
  double alpha_s;
  double alpha_min;
  double r_max;
  double lambda_bym;

  //******************************************************************************************
  //  EOS Config parameters
  //******************************************************************************************
  double a_trento;
  double s_chop;
  vector<SplineSet> eos_interped;

  //******************************************************************************************
  //  Multiple Use Config parameters
  //******************************************************************************************
  double grid_max;
  double grid_step;
  double tau_0; //  EOS, Event, Splitting
  double e_thresh; // Splitting and event
  string charge_type; //  Flag for tracking BSQ or UDS charges

  int tracked_charge;
  int current_event;  //  Tracks what event is being handled
  int grid_points;  // This is calculated to be 2*(grid_max/grid_step)
  //#CONFIGPARAM
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Internal Functions
//##########################################################################################
  //  Initializes map used to read config file
  void Initialize();

  //  Copy function for IO class, called by operator= and implicit copy functions
	void CopyIO(const IO &e);

  // Convert event input to energy
  void ConvertEvent(vector<vector<double>> &input, double &total);

  //  Definitions of possible output formats
  void OutputConfig(string file_name);  //  Prints copy of config file to output directory for reference
  void OutputFullDensityGrids(vector<vector<double>> &density_grid, string file_name); //  Prints full density grids with 0s
  void OutputSparseDensityGrids(vector<vector<double>> &density_grid, string file_name);  //  Only prints valued points
  void OutputSparseDensityGrids(vector<vector<vector<double>>> &density_grid, double tot_energy, string file_name);  //  Only prints valued points
  void OutputEccentricities(double total_entropy, vector<vector<double>> eccentricities, string density_type, string file_name); //  Prints eccentricities
  void OutputQuarkCounts(double total_entropy, int gluon, int up, int down, int strange, int charm, string file_name); //  Prints quark counts
//__________________________________________________________________________________________

public:

//__________________________________________________________________________________________
//##########################################################################################
//  Basic Class Functions
//##########################################################################################
  IO(string configFile);  //  Class constructor, must specify path to configuration file
  ~IO();  //  Destructor, used to clear all data stored by class

  IO(const IO &original); //  Implicit copy function, newIOObject(oldIOObject)
  IO& operator=(const IO& original);  //  Defines what happens when you use = operator on class
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  IO Specific Functions
//##########################################################################################
  // Initialize the Event Object
  Event InitializeEvent();
  // Initialize the Splitter Object
  Splitter InitializeSplitter();
  // Initialize equation of state
  void InitializeEOS();

  string GetOutputDir() { return output_dir;  }
  string GetTest() { return test_;  }

  //  Read single event
  Event ReadEvent(Event event_in);

  //  Writes single event to file
  void WriteEvent(Event event);

  //  Test for end of event list
  bool LastEvent();
//__________________________________________________________________________________________
};
#endif
