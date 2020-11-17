#include "io.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Reads in config file and initializes variables
//##########################################################################################
IO::IO(string configFile)
{
  Initialize(); // Initialze config params and map for reading in config params

  //  config input stream
  ifstream input;
  input.open(configFile);

  string var_type;  //  Used to record parameter type read from file
  string data;
  string variabletype;
  //******************************************************************************************
  //  Loop through config file
  //******************************************************************************************
  while (!input.eof())
	{
    input >> var_type;  //  Read parameter type

    //  Switch through the possible parameter types
    //  uses var_type as key to map then reads value to class variable
    switch(mapConfigParams[var_type])
    {
      //  case ConfigParam (Does var_type map to ConfigParam?)
      //  input >> input_var (Read in value of var_type)
      //  break; (Stop checking switch and move on)
      case eventinputfile: input >> event_input_file;
      case outputfolder: input >> output_folder;

      case dataformat:
        getline(input, data);
        istringstream datatypes(data);
        int num = 0;
        while (datatypes >> variabletype)
        {
          if (variabletype == "ev") data_locations[0] = num;
          if (variabletype == "b") data_locations[1] = num;
          if (variabletype == "Npart") data_locations[2] = num;
          if (variabletype == "Mult") data_locations[3] = num;
          if (variabletype == "s") data_locations[4] = num;
          if (variabletype == "e2") data_locations[5] = num;
          if (variabletype == "phi2") data_locations[6] = num;
          if (variabletype == "e3") data_locations[7] = num;
          if (variabletype == "phi3") data_locations[8] = num;
          if (variabletype == "e4") data_locations[9] = num;
          if (variabletype == "phi4") data_locations[10] = num;
          if (variabletype == "e5") data_locations[11] = num;
          if (variabletype == "phi5") data_locations[12] = num;
          if (variabletype == "rad") data_locations[13] = num;
          num++;
        }
      //#CONFIGPARAM
    }// End of switch
  }// End of while loop
cout << output_folder << endl;
  input.close();
//  OutputConfig(output_dir + "run_parameters" + to_string(current_event) + ".dat");
}// End of Class constructor
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class Destructor
//##########################################################################################
IO::~IO()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  IO Copy Function
//##########################################################################################
void IO::CopyIO(const IO &e)
{
  event_input_file = e.event_input_file;
  output_folder = e.output_folder;
  data_locations = e.data_locations;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
IO::IO(const IO &original)
{
  CopyIO(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
IO& IO::operator= (const IO& original)
{
	CopyIO(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Initialze variables and map
//##########################################################################################
void IO::Initialize()
{
  //******************************************************************************************
  //  Set variables to default values
  //******************************************************************************************
  event_input_file = "";
  output_folder = "";
  data_locations.resize(14, -1);
  //******************************************************************************************
  //  Initialze map for reading in config file
  //******************************************************************************************
  mapConfigParams["event_input_file"] = eventinputfile;
  mapConfigParams["output_folder"] = outputfolder;
  mapConfigParams["data_format"] = dataformat;
  //#CONFIGPARAM
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Copy config file to ICCING output folder for future reference
//##########################################################################################
void IO::OutputConfig(string file_name)
{
  ofstream output;
  output.open(file_name);

    //#CONFIGPARAM

    output.close();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
// Read event data
//##########################################################################################
vector<Event> IO::ReadEvents()
{
  //  Input file stream
  ifstream input;
  input.open(event_input_file);
  cout << event_input_file << endl;

  string input_data;
  vector<double> event_data;
  int x;
  event_data.resize(data_locations.size(), 0.0);

  vector<Event> event_list;

  while(getline(input, input_data))
  {
    istringstream split_data(input_data);
    x = 0;
    while (split_data >> event_data[x]) { x++;  }

    Event current_event;
    if (data_locations[0] != -1) current_event.event_num = event_data[data_locations[0]];
    if (data_locations[1] != -1) current_event.impact_parameter = event_data[data_locations[1]];
    if (data_locations[2] != -1) current_event.number_of_participants = event_data[data_locations[2]];
    if (data_locations[3] != -1) current_event.multiplicity = event_data[data_locations[3]];
    if (data_locations[4] != -1) current_event.entropy = event_data[data_locations[4]];
    if (data_locations[5] != -1) current_event.eccentricity[2] = event_data[data_locations[5]];
    if (data_locations[6] != -1) current_event.phi[2] = event_data[data_locations[6]];
    if (data_locations[7] != -1) current_event.eccentricity[3] = event_data[data_locations[7]];
    if (data_locations[8] != -1) current_event.phi[3] = event_data[data_locations[8]];
    if (data_locations[9] != -1) current_event.eccentricity[4] = event_data[data_locations[9]];
    if (data_locations[10] != -1) current_event.phi[4] = event_data[data_locations[10]];
    if (data_locations[11] != -1) current_event.eccentricity[5] = event_data[data_locations[11]];
    if (data_locations[12] != -1) current_event.phi[5] = event_data[data_locations[12]];
    if (data_locations[13] != -1) current_event.radius = event_data[data_locations[13]];

    event_list.push_back(current_event);
  }

  input.close();

  return event_list;
}
//__________________________________________________________________________________________

void IO::OutputObservables(vector<vector<double>> observables, string file)
{
  ofstream output;
  output.open(output_folder + file);
cout << output_folder + file << endl;
  for (int i = 0; i < observables.size(); i++)
  {
    for (int j = 0; j < observables[i].size(); j++)
    {
      output << observables[i][j];
    }
    output << endl;
  }
  output.close();
}
