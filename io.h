#ifndef IO_H
#define IO_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################
#include "event.h"
#include "ecc.h"
#include "splitting.h"
#include "probdist.h"
//__________________________________________________________________________________________

using namespace std;

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

    eosemmitlines,
    eosscol,
    eosecol,
    atrento,
    echop,

    bmin,
    bmax,
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
  vector<vector<double>> flavor_chemistry;  //  flavor_chemistry[Qs(GeV)][quark_prob] quark_prob: 0 = Qs, 1 = up, 2 = down ...
  string dipole_model;  //  Future Update
  double alpha_s;
  double alpha_min;
  double r_max;

  //******************************************************************************************
  //  EOS Config parameters
  //******************************************************************************************
  int eos_emmit_lines;
  int eos_s_col;
  int eos_e_col;
  double a_trento;
  double e_chop;


  //******************************************************************************************
  //  Multiple Use Config parameters
  //******************************************************************************************
  int b_min;
  int b_max;
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

  //  Definitions of possible output formats
  void OutputConfig(string file_name);  //  Prints copy of config file to output directory for reference
  void OutputFullDensityGrids(vector<vector<double>> density_grid, string file_name); //  Prints full density grids with 0s
  void OutputSparseDensityGrids(vector<vector<double>> density_grid, string file_name);  //  Only prints valued points
  void OutputEccentricities(Eccentricity &ecc, string file_name); //  Prints eccentricities
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
  Event InitializeEvent();
  Splitter InitializeSplitter();
  ProbDist InitializeProbDist();
//  EOS InitializeEOS();

  Event ReadEvent();  //  Read single event

  void WriteEvent(Event event); //  Writes single event to file

  bool LastEvent(); //  Test for end of event list
//__________________________________________________________________________________________
};
#endif
