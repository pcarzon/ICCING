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
#include <sstream>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################
#include "event.h"
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
    eventinputfile,
    outputfolder,
    dataformat
    //#CONFIGPARAM
  };

  //  Maps values of ConfigParams to keyword strings used in config file
  map<string, ConfigParams> mapConfigParams;

  //  File variables
  string event_input_file;
  string output_folder;

  //******************************************************************************************
  //  Internal Config parameters
  //******************************************************************************************
  vector<int> data_locations;
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
vector<Event> ReadEvents();

void MakeDirectory(int bin);

void OutputObservables(vector<vector<double>> observables, string file);
//__________________________________________________________________________________________
};
#endif
