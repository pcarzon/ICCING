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
      case dataformat:
        getline(input, data);
        istringstream datatypes(data);
        int num = 0;
        while (datatypes >> variabletype)
        {
          if (variabletype == "ev") { cout << "got event number" << endl; data[0] = num;}
          if (variabletype == "s") data[1] = num;
          if (variabletype == "e2") data[2] = num;
          if (variabletype == "phi2") data[3] = num;
          if (variabletype == "e3") data[4] = num;
          if (variabletype == "phi3") data[5] = num;
          if (variabletype == "e4") data[6] = num;
          if (variabletype == "phi4") data[7] = num;
          if (variabletype == "e5") data[8] = num;
          if (variabletype == "phi5") data[9] = num;
          if (variabletype == "rad") data[10] = num;
          num++;
        }
        cout << "data size " << data.size() << endl;

        for (int i = 0; i < data.size(); i++)
        cout << i << " " << data[i] << endl;
        exit(0);
      //#CONFIGPARAM
    }// End of switch
  }// End of while loop

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
  data = e.data;
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
  data.resize(11, 0);
  cout << "data size " << data.size() << endl;
  //******************************************************************************************
  //  Initialze map for reading in config file
  //******************************************************************************************
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
