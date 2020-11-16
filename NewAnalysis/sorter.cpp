#include "sorter.h"

Sorter::Sorter(vector<Event> all_events)
{
  sort(all_events.begin(), all_events.end(), compareByEntropy);
  sorted_events = all_events;

  events_1percent_centrality_bin = int(sorted_events.size()/100);
}

/*double Sorter::GenerateCentralityCuts(string folder, string sort_file, vector<Event> all_events)
{

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax, cmin;


  sort(all_events.begin(), all_events.end(), compareByEntropy);

  int events_per_bin = int (all_events.size()/100);
  //cout << events_per_bin << " " << all.size() << endl;

	vector<double> cuts;
	for (int c = 0; c < 100; c++)
	{
    cuts.push_back(all_events[c*events_per_bin].entropy);
  }

	int end = events_per_bin*100 - 1;

	cuts[100] = all_events[end].entropy;


  string name = folder + "/" + sort_file + ".dat";
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name << endl;
 		exit(1);
 	}

  for (int c = 0; c < 99; c++)
	{
		PRINT << c << " " << c+1 << " " <<  cuts[c]  << " "  << cuts[c + 1] << endl;
	}

  PRINT.close();
  double cuts1 = cuts[1];


  string name2 = folder + "/" + sort_file + "5.dat";
  ofstream PRINT2(name2.c_str());
	if (!PRINT2.is_open())
 	{
 		cout << "Can't open " << name2 << endl;
 		exit(1);
 	}

  for (int c = 0; c < 99; c = c + 5)
	{
		PRINT2 << c << " " << c + 5 << " " <<  cuts[c]  << " "  << cuts[c + 5] << endl;
	}

  PRINT2.close();

  double cuts5 = cuts[5];

  string name3 = folder + "/" + sort_file + "10.dat";
  ofstream PRINT3(name3.c_str());
	if (!PRINT3.is_open())
 	{
 		cout << "Can't open " << name3 << endl;
 		exit(1);
 	}

  for (int c = 0; c < 99; c = c + 10)
	{
		PRINT3 << c << " " << c + 10 << " " <<  cuts[c]  << " "  << cuts[c + 10] << endl;
	}

  double cuts10 = cuts[10];
  PRINT3.close();

  //cout << cuts[1] << " " << cuts10 << " " << cuts1 << endl;

  return cuts10;

}

int Sorter::GenerateNumberOfParticipantsCuts(string folder,vector<Event> all_events)
{

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax,cmin;


  sort (all_events.begin(), all_events.end(), compareByNpart);

  int events_per_bin = all_events.size()/100;

	vector<int> cuts;
	for (int c = 0; c < 100; c++)
	{
		cuts.push_back(all_events[c*events_per_bin].number_of_participants);
	}

	int end = all_events.size() - 1;
	cuts[100] = all_events[end].number_of_participants;



  string name = folder + "/npartcen1.dat";
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name << endl;
 		exit(1);
 	}

  for (int c = 0; c < 99; c++)
	{
		PRINT << c << " " << c + 1 << " " <<  cuts[c]  << " "  << cuts[c + 1] << endl;
	}
  double cuts1 = cuts[1];
  PRINT.close();


  string name2 = folder + "/npartcen5.dat";
  ofstream PRINT2(name2.c_str());
	if (!PRINT2.is_open())
 	{
 		cout << "Can't open " << name2 << endl;
 		exit(1);
 	}

  for (int c = 0; c < 99;c = c + 5)
	{
		PRINT2 << c << " " << c + 5 << " " <<  cuts[c]  << " "  << cuts[c + 5] << endl;
	}
  double cuts5 = cuts[5];
  PRINT2.close();

  string name3 = folder + "/npartcen10.dat";
  ofstream PRINT3(name3.c_str());
	if (!PRINT3.is_open())
 	{
 		cout << "Can't open " << name3 << endl;
 		exit(1);
 	}

  for (int c = 0; c < 99; c = c + 10)
	{
		PRINT3 << c << " " << c + 10 << " " <<  cuts[c]  << " "  << cuts[c + 10] << endl;
	}
  double cuts10 = cuts[10];
  PRINT3.close();

	//cout << cuts[1] << " " << cuts10 << " " << cuts1 << endl;

  return cuts10;

}
*/
//__________________________________________________________________________________________
//##########################################################################################
//  Class Destructor
//##########################################################################################
Sorter::~Sorter()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Sorter Copy Function
//##########################################################################################
void Sorter::CopySorter(const Sorter &e)
{
  sorted_events = e.sorted_events;
  events_1percent_centrality_bin = e.events_1percent_centrality_bin;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Sorter::Sorter(const Sorter &original)
{
  CopySorter(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Sorter& Sorter::operator= (const Sorter& original)
{
	CopySorter(original);
	return *this;
}
//__________________________________________________________________________________________

vector<vector<Event>> Sorter::SortEccentricitiesIntoCentralityBins(int bin_width)
{
  vector<vector<Event>> binned_events;
  binned_events.resize(100/bin_width);

  for (int bin = 0; bin < binned_events.size(); bin++)
  {
    double average_entropy = 0;

    for (int i = 0; i < events_1percent_centrality_bin*bin_width; i++)
    {
      binned_events[bin].push_back(sorted_events[events_1percent_centrality_bin*bin + i]);
      average_entropy += sorted_events[events_1percent_centrality_bin*bin + i].entropy;
    }
    cout << "bin " << bin << " entropy " << average_entropy/binned_events[bin].size() << " lower bound " << binned_events[bin][binned_events.size()].entropy << endl;
  }
//  cout << "events per bin " << events_1percent_centrality_bin*bin_width << endl;
//  for (int bin = 0; bin < binned_events.size(); bin++)
//  {
  //  cout << "bin " << bin << " size = " << binned_events[bin].size() << endl;
//  }
  return binned_events;
}
