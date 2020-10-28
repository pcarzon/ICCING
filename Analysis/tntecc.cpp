.entropysorted_eccentricities#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string.h>
#include <vector>
#include <algorithm>
#include <tuple>
#include <iterator>



//program specific classes
#include "sc.h"

//namespace
using namespace std;

//*******************************************************************************************************************************************************************************
//Declaration of functions
//*******************************************************************************************************************************************************************************

// sort - sorting events into centrality bins (actually does the work) (sorted_eccentricities[centrality][i])
void SortEccentricitiesIntoCentralityBins(string name,string sfile,vector<eccs> all_eccentricities,  vector< vector<eccs> > & out, vector<double> & cens,int bin);

// NPsc - symetric cummulants without fancy stuff as one does in mathematica
void NPsc(string name2,  string chg,pars p,  vector<eccs>  & ec);
//normalized sc

// bsort - impact parameter
void bsort(string folder, vector<eccs> all_eccentricities);

// MaximumNumberOfParticipants - Checks all events for the one with the greatest number of participants
int MaximumNumberOfParticipants( vector<eccs> s1)
{
	int ss=s1.size();
	int max=0;
	for (int i=0;i<ss;i++)
	{
		if (s1[i].number_of_participants>max) max=s1[i].number_of_participants;
	}
  return max;
}

// pred - predictions (with hydro run)
double pred(double k1, double k2, double ec);

// PredictCummulants - predict v22 and v24 using kappas (skanda's work) comment out for now
void PredictCummulants(string ictype,string runtype, string runname,string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens, int bin);

// GenerateCentralityCuts - sort in centrality (make file of cuttoffs)
double GenerateCentralityCuts(string folder, string sort_file, vector<eccs> all_eccentricities);

// GenerateNumberOfParticipantsCuts - sort in npart (make file of cuttoffs)
int GenerateNumberOfParticipantsCuts(string folder, vector<eccs> all_eccentricities);

//void sort2(string name,vector<eccs> all,  vector< vector<eccs> > & out, vector<double> & cens);

//void sort3(string name,vector<eccs> all,  vector< vector<eccs> > & out, vector<double> & cens);

// ultracen -
void ultracen(string name2, string chg, vector<eccs> all_eccentricities, double centrality_cuts, int number_of_participants_cuts);

// ucprintCGC - uses CGC from trento
void ucprintCGC(string name2,string chg,string pname,vector<eccs> ssort);

// ucprint -
void ucprint(string name2,string chg,string pname,vector<eccs> ssort);

// vCGC - CGC from trento accross all centralities
void vCGC(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens,int bin);

// ReadEccentricities - read in eccentricities from folder (add an if to skip nans)
string ReadEccentricities(string name2, vector<eccs> & all_eccentricities);

// vns - calculate all Vnm's with error bars (easy mode with statistics)
void CalculateCummulants(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens,int bin);

// vnsnorat - same as vns but with out stuff like v24/v22
void CalculateCummulantsWithoutRatios(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens,int bin);
//void vnsb1(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens);

//void vns2(string name2,string chg,  vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens);

//void vns3(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens);

// radii - average radius vs centrality
void radii(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens,int bin);


// compareByLength - compare by multiplicity
bool compareByLength(const eccs &a, const eccs &b)
{
    return a.entropy > b.si;
}

// compareByNpart -
bool compareByNpart(const eccs &a, const eccs &b)
{
    return a.number_of_participants > b.number_of_participants;
}
//*******************************************************************************************************************************************************************************
//End of function declarations
//*******************************************************************************************************************************************************************************

//*******************************************************************************************************************************************************************************
//main
//*******************************************************************************************************************************************************************************
int main (int argc, char *argv[]) //execute with ./a.out folder list_name sortfile_name PID
{

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// folder_path = folder path to file of eccentricities (assumption /trento/folder_path)
	// sname = name of list of eccentricities without .dat extension
	// chg = set to 0
	// sfile = output file name
	// runtype = folder for running with kappas (from Skanda)
	// runname = kappa file
	string folder_path, sname, chg, sfile, runtype, runname;

	if (argv[1])	{	folder_path = argv[1];	}
	else	{	folder_path = "";	}

	if (argv[2])	{	sname = argv[2];	}
	else	{	sname = "";	}

	if (argv[3])	{	sfile = argv[3];	}
	else	{	sfile = "";	}

	if (argv[4])	{	chg = argv[4];	}
	else	{	chg = "0";	}

	//specify IC folder to obtain kappas
  if (argv[5])	{	runtype = argv[5];	}
	else	{	runtype = "0";	}

	//specify IC folder to obtain kappas
  if (argv[6])	{	runname = argv[6];	}
	else	{	runname = "0";	}
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//	Read Eccentricities from file
	vector<eccs> all_eccentricities;
	string input_eccentricities_folder = folder_path + "/" + sname;
	string data_type = ReadEccentricities(input_eccentricities_folder, all_eccentricities);
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//	Generate Cuts for sorting eccentricities
  double centrality_cuts = GenerateCentralityCuts(folder_path, sfile, all_eccentricities);
  int number_of_participants_cuts = GenerateNumberOfParticipantsCuts(folder_path, all_eccentricities);
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	cout << all_eccentricities.size() << endl;
	vector< vector<eccs> > sorted_eccentricities;


	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//	Sort and Calculate in 5% centrality bins
	vector<double> CentralityBins_5Percent;
  bsort(folder_path, all_eccentricities);
	SortEccentricitiesIntoCentralityBins(folder_path, sfile, all_eccentricities, sorted_eccentricities, CentralityBins_5Percent, 5); //5% bins
	int ssize = sorted_eccentricities.size();
	//	for (int j=0;j<ssize;j++ ) SC.correct(sorted_eccentricities[j]);

	CalculateCummulants(folder_path, chg, sorted_eccentricities, CentralityBins_5Percent, 5);
	if (data_type == "Trento + CGC + R")
	{	radii(folder_path, chg, sorted_eccentricities, CentralityBins_5Percent, 5);	}
	CalculateCummulantsWithoutRatios(folder_path, chg, sorted_eccentricities, CentralityBins_5Percent, 5);
	PredictCummulants(folder_path, runtype, runname, chg, sorted_eccentricities, CentralityBins_5Percent, 5);
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//	Sort and Calculate in 1% centrality bins
	vector<double> CentralityBins_1Percent;
	sorted_eccentricities.clear();
	SortEccentricitiesIntoCentralityBins(folder_path, sfile, all_eccentricities, sorted_eccentricities, CentralityBins_1Percent, 1);
  if (data_type == "Trento + CGC" || data_type == "Trento + CGC + R")
	{	vCGC(folder_path, chg, sorted_eccentricities, CentralityBins_1Percent, 1);	}

  CalculateCummulants(folder_path, chg, sorted_eccentricities, CentralityBins_1Percent, 1);
  if (data_type == "Trento + CGC + R")
	{	radii(folder_path, chg, sorted_eccentricities, CentralityBins_1Percent, 1);	}
  CalculateCummulantsWithoutRatios(folder_path, chg, sorted_eccentricities, CentralityBins_1Percent, 1);
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//	Calculate UltraCentral observables
  if (data_type == "Trento + CGC" || data_type == "Trento + CGC + R")
	{	ultracen(folder_path, chg, all_eccentricities, centrality_cuts, number_of_participants_cuts);	}
	//  PredictCummulants(folder_path,runtype,runname,chg,sorted_eccentricities,CentralityBins_1Percent,1);
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Calculate Symmetric Cummulants
	sc SC;
	string pre = folder_path + "/obsb3_" + chg + "/";
	string into;

	pars p1;
	p1.n = 3;
	p1.m = 2;
  NPsc(folder_path, chg, p1, all_eccentricities);
	into = pre + "SCM32comb";
	SC.runec(into, p1, sorted_eccentricities);

	p1.n = 4;
	p1.m = 2;
	into = pre + "SCM42comb";
	SC.runec(into, p1, sorted_eccentricities);

	p1.n = 4;
	p1.m = 3;
	into = pre + "SCM43comb";
	SC.runec(into, p1, sorted_eccentricities);
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
}
//*******************************************************************************************************************************************************************************
//End of main
//*******************************************************************************************************************************************************************************


//*******************************************************************************************************************************************************************************
//Function Definitions
//*******************************************************************************************************************************************************************************
string ReadEccentricities(string input_eccentricities_folder, vector<eccs> & all_eccentricities)
{
	int num_data_entrys;
	string data_type = "Default Trento";

	string input_eccentricities_file = input_eccentricities_folder + ".dat";
	ifstream input(input_eccentricities_file.c_str());
	if (!input.is_open())
 	{
 		cout << "Can't open " << input_eccentricities_file << endl;
 		exit(1);
 	}

	string line;
	while (getline(input, line))
	{

		vector<double> y (15, 0) ;
		vector<string> x = split(line, ' ');

    int tsize = x.size();
  	if (tsize < 8) break;

    num_data_entrys = 8;

    eccs current_event;
		for(int j = 0; j < num_data_entrys; j++)
		{

			//	y[j] = stod(x[j]);
			stringstream s;
			s << x[j];
			s >> y[j];
			//	cout << y[j] << " " << endl;
		}

		if (num_data_entrys == 5)
		{// old mckln files
			current_event.entropy= y[0];
			current_event.multiplicity = y[0];
	  	current_event.number_of_participants = y[0];
			current_event.eccentricity[2] = y[1];
			current_event.eccentricity[3] = y[3];
			current_event.eccentricity[4] = 0;
			current_event.eccentricity[5] = 0;
			current_event.flow_harmonic[2] = y[1];
			current_event.flow_harmonic[3] = y[3];
			current_event.flow_harmonic[4] = 0;
    	current_event.impact_parameter = 0;
			data_type = "Old MCKLN";
		}

		if (num_data_entrys > 11)
		{//original trento output
 			current_event.impact_parameter = y[1];
			current_event.entropy= y[3];
			current_event.multiplicity = y[3];
 			current_event.number_of_participants = y[2];
			current_event.eccentricity[2] = y[4];
			current_event.eccentricity[3] = y[5];
			current_event.eccentricity[4] = y[6];
			current_event.eccentricity[5] = y[7];
			current_event.flow_harmonic[2] = y[4];
			current_event.flow_harmonic[3] = y[5];
			current_event.flow_harmonic[4] = y[6];
			data_type = "ICCING Output";
		}
		//	cout << "3" << endl;
		//work here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
		if (num_data_entrys > 5)
		{//original trento output
 			current_event.impact_parameter = y[1];
			current_event.entropy= y[3];
			current_event.multiplicity = y[3];
 			current_event.number_of_participants = y[2];
			current_event.eccentricity[2] = y[4];
			current_event.eccentricity[3] = y[5];
			current_event.eccentricity[4] = y[6];
			current_event.eccentricity[5] = y[7];
			current_event.flow_harmonic[2] = y[4];
			current_event.flow_harmonic[3] = y[5];
			current_event.flow_harmonic[4] = y[6];
		}

		//this is all of trento stuff
  	if (num_data_entrys > 8)
		{// trento+CGCoutput
			current_event.flow_harmonics_no_decays[1] = y[8];//equivalent to "multiplicity"
  		current_event.flow_harmonics_no_decays[2] = pow(y[8], 1)*y[9];
			current_event.flow_harmonics_no_decays[3] = pow(y[8], 1)*y[10];
			current_event.flow_harmonics_no_decays[4] = pow(y[8], 1)*y[11];
			data_type = "Trento + CGC";
			if (num_data_entrys > 13)
	 		{
		 		current_event.radius = y[13];
	   		data_type = "Trento + CGC + R";
 	 		}
		}
		//cout << "try" << endl;
		all_eccentricities.push_back(current_event);
	}

	input.close();

	return data_type;
}

void bsort(string folder, vector<eccs> all_eccentricities){

  double becc2[12] = {0}, becc3[12] = {0}, cnt[12] = {0};

  int tot = all_eccentricities.size();
	for (int i = 0; i < tot; i++)
	{
		for (int b = 0; b < 12; b++)
		{
      int bmin = b - 1;
      if (bmin < 0) bmin = 0;
			if (all_eccentricities[i].impact_parameter <= b && all_eccentricities[i].impact_parameter >= bmin)
			{
        becc2[b] += all_eccentricities[i].eccentricity[2];
        becc3[b] += all_eccentricities[i].eccentricity[3];
        cnt[b]++;
        break;
      }
		}
	}

  for (int b = 0; b < 12; b++)
	{
    cout << b << " " << cnt[b] << " " << becc2[b]/cnt[b] << " " << becc3[b]/cnt[b] << endl;
  }
}


double GenerateCentralityCuts(string folder, string sort_file, vector<eccs> all_eccentricities)
{

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax, cmin;


	cout << "start sort" << endl;
  sort (all_eccentricities.begin(), all_eccentricities.end(), compareByLength);

  int dbin = int (all_eccentricities.size()/100);
  //cout << dbin << " " << all.size() << endl;

	vector<double> cuts;
	for (int c = 0; c < 100; c++)
	{
    cuts.push_back(all_eccentricities[c*dbin].si);
  }

	int end = dbin*100 - 1;

	cuts[100] = all_eccentricities[end].si;


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

int GenerateNumberOfParticipantsCuts(string folder,vector<eccs> all_eccentricities)
{

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax,cmin;


  sort (all_eccentricities.begin(), all_eccentricities.end(), compareByNpart);

  int dbin = all_eccentricities.size()/100;

	vector<int> cuts;
	for (int c = 0; c < 100; c++)
	{
		cuts.push_back(all_eccentricities[c*dbin].number_of_participants);
	}

	int end = all_eccentricities.size() - 1;
	cuts[100] = all_eccentricities[end].number_of_participants;



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

void SortEccentricitiesIntoCentralityBins(string folder, string sort_file, vector<eccs> all_eccentricities,  vector< vector<eccs> > & out, vector<double> & cens,int bin)
{

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax, cmin;


	for (int i = 0; i < (100 - bin); i += bin)
	{
		cmax.push_back(i + bin);
		cmin.push_back(i);
		cens.push_back(i + bin/2.);
	}

	cout << "check sizes cmin " << cmin.size() << " " << cmax.size() << endl;

  int censize = cens.size();
	vector<double> nph, npl;
	nph.resize(censize);
	npl.resize(censize);


	string named = folder + "/" + sort_file + ".dat";
	//cout << "named  " << named << endl;
	ifstream inputd(named.c_str());
	if (!inputd.is_open())
 	{
 		cout << "Can't open " << named << endl;
 		exit(1);
 	}
 	string line;
 	int imin = 0,imax = 0;

	while (getline(inputd, line))
	{
		vector<double> y (2, 0) ;
		vector<double> y2 (4, 0) ;
		vector<string> x = split(line, ' ');
		if (x.size()<4) vector<string> x = split(line, '\t');

		//cout << "x.size  " << x.size() << endl;

		for(int j = 0; j < 2; j++)
		{
			stringstream s;
			s << x[j];
			s >> y[j];
		}
		for(int j = 2; j < 4; j++)
		{
			stringstream s;
			s << x[j];
			s >> y2[j];
		}

		//cout <<"imin  " << imin << " imax  " << imax  << " nph  " << nph.size() << " npl  " << npl.size() << endl;
		if (cmin[imin] == y[0])
		{//min of the centrality bin
	 		nph[imin] = y2[2]; //max entropy for centrality bin
	 		imin++;
	 		if (imin > censize) break;
	 	}
		if (cmax[imax] == y[1])
		{ //max of the centrality bin
	 		npl[imax] = y2[3];//min entropy for centrality bin
	 		imax++;
	 		if (imax > censize) break;
	 	}

	}
	inputd.close();


	//sets up variables to the approparite number of bins in centrality
	int tot = all_eccentricities.size();

	out.resize(censize);
	for (int i = 0; i < tot; i++)
	{
		for (int c = 0; c < censize; c++)
		{
			if (all_eccentricities[i].entropy<= nph[c] && all_eccentricities[i].entropy> npl[c])
			{
				out[c].push_back(all_eccentricities[i]);
				break;
			}
		}
	}

}


void NPsc(string name2,  string chg,pars p,  vector<eccs>  & ec)
{
  // sort eccentricities by Npart, store in list

  int max = MaximumNumberOfParticipants(ec);
  int size = ec.size();

  vector<vector<eccs>> list;
  list.resize(max + 1);
  for (int i = 0; i < size; i++)
	{
  	list[ec[i].number_of_participants].push_back(ec[i]);
  }

  vector<double> sc32, sc42, sc43, v2, v3, v4;
  sc32.resize(max + 1, 0);
  sc42.resize(max + 1, 0);
  sc43.resize(max + 1, 0);
  v2.resize(max + 1, 0);
  v3.resize(max + 1, 0);
  v4.resize(max + 1, 0);

  for (int n = 0; n <= max; n++)
	{
  	int full = list[n].size();
  	for (int f = 0; f < full; f++)
		{
  		v2[n] += pow(list[n][f].eccentricity[2], 2);
  		v3[n] += pow(list[n][f].eccentricity[3], 2);
  		v4[n] += pow(list[n][f].eccentricity[4], 2);
  		sc32[n] += pow(list[n][f].eccentricity[2]*list[n][f].eccentricity[3], 2);
  		sc42[n] += pow(list[n][f].eccentricity[2]*list[n][f].eccentricity[4], 2);
 			sc43[n] += pow(list[n][f].eccentricity[4]*list[n][f].eccentricity[3], 2);
  	}
  	v2[n] /= full;
  	v3[n] /= full;
  	v4[n] /= full;
  	sc32[n] /= full;
  	sc42[n] /= full;
  	sc43[n] /= full;
  }
  //cout << "list size 0 " << list[0].size() << endl << "list size 1 " << list[1].size() << endl << "list size max " << list[max].size() << endl << "list size " << list.size() << endl;
  string name = name2 + "/obsb3_" + chg + "/eNSC_Npart.dat";
  cout << "max " << max << endl;
  ofstream fileout(name.c_str());
  cout << "Print is open" << endl;
	if (!fileout.is_open())
 	{
 		cout << "Can't open " << name << endl;
 		exit(1);
 	}
	//cout << "testing" << endl;
	for (int n = 1; n <= max; n++)
	{
		//cout << n << endl;
    fileout << n << " " << sc32[n]/(v2[n]*v3[n]) - 1  << " " << sc42[n]/(v2[n]*v4[n]) - 1  << " "<< sc43[n]/(v4[n]*v3[n]) - 1  <<  endl;
  }

  fileout.close();
	cout << "done np sc" << endl;
}


void ultracen(string name2,  string chg,vector<eccs> all_eccentricities, double centrality_cuts, int number_of_participants_cuts)
{
  int tot = all_eccentricities.size();
  vector<eccs> ssort, npsort;
  for (int i = 0; i < tot; i++)
	{
    if (all_eccentricities[i].si >= centrality_cuts) ssort.push_back(all_eccentricities[i]);
    if (all_eccentricities[i].number_of_participants >= number_of_participants_cuts) npsort.push_back(all_eccentricities[i]);
  }

	cout << ssort.size() << " " << npsort.size() << endl;

  SortEccentricitiesIntoCentralityBins(ssort.begin(), ssort.end(), compareByLength);
  SortEccentricitiesIntoCentralityBins(npsort.begin(), npsort.end(), compareByLength);

 	string sn = "entropy";
  ucprint(name2, chg, sn, ssort);
  sn = "npart";
  ucprint(name2, chg, sn, npsort);
  sn = "entropy";
  ucprintCGC(name2, chg, sn, ssort);
  sn = "npart";
  ucprintCGC(name2, chg, sn, npsort);
	cout << "mpart" << endl;
}

void ucprint(string name2, string chg, string pname, vector<eccs> ssort)
{
  int slen = ssort.size();

  int sbin = slen/20;

  int sc = 0;
  int npc = 0;
  cout << "starting" << endl;
  vector<double> v2, v3, sc32, M;
  double mv2 = 0, mv3 = 0, msc32 = 0, mM = 0;
  while(sc < slen)
	{
    double c2 = 0, c3 = 0, c32 = 0, cM = 0;
    for(int i = sc; i < (sc + sbin); i++)
		{
      if (i >= slen) break;
      c2 += pow(ssort[i].flow_harmonics[2], 2);
      c3 += pow(ssort[i].flow_harmonics[3], 2);
      c32 += pow(ssort[i].flow_harmonics[2]*ssort[i].flow_harmonics[3], 2);
      cM += ssort[i].multiplicity;
    }
    mv2 += c2;
    mv3 += c3;
    msc32 += c32;
    mM += cM;
    v2.push_back(sqrt(c2/sbin));
    v3.push_back(sqrt(c3/sbin));
    sc32.push_back(c32/sbin);
    M.push_back(cM/sbin);

  	sc += sbin;
  }

  double c2 = 0, c3 = 0, c32 = 0, cM = 0;
  for(int i = sc; i < slen; i++)
	{
    c2 += pow(ssort[i].flow_harmonics[2], 2);
    c3 += pow(ssort[i].flow_harmonics[3], 2);
    c32 += pow(ssort[i].flow_harmonics[2]*ssort[i].flow_harmonics[3], 2);
    cM += ssort[i].multiplicity;
  }
  mv2 += c2;
  mv3 += c3;
  msc32 += c32;
  mM += cM;
  v2.push_back(sqrt(c2/sbin));
  v3.push_back(sqrt(c3/sbin));
  sc32.push_back(c32/sbin);
  M.push_back(cM/sbin);


  mv2 /= slen;
  mv3 /= slen;
  msc32 /= slen*mv2*mv3;
  mM /= slen;

  cout << slen << endl;

  int olen = v2.size();
  vector<double> erv2, erv3, ersc32;
  for(int i = 0; i < olen; i++)
	{
    double sv2 = 0, sv3 = 0, ssc32 = 0;


    for(int j = i*sbin; j < ((i + 1)*sbin); j++)
		{
      if (j >= slen) break;
      sv2 += pow(v2[i]/sqrt(mv2) - sqrt((pow(v2[i], 2)*sbin - pow(ssort[j].flow_harmonics[2], 2))/(sbin - 1))/sqrt(mv2), 2);
      sv3 += pow(v3[i]/sqrt(mv3) - sqrt((pow(v3[i], 2)*sbin - pow(ssort[j].flow_harmonics[3], 2))/(sbin - 1))/sqrt(mv3), 2);
      ssc32 += pow(sc32[i]/(pow(v2[i], 2)*pow(v3[i], 2)) - (sbin - 1)*(sc32[i]*sbin - pow(ssort[i].flow_harmonics[2]*ssort[i].flow_harmonics[3], 2))/((pow(v2[i], 2)*sbin - pow(ssort[j].flow_harmonics[2], 2))*(pow(v3[i], 2)*sbin - pow(ssort[j].flow_harmonic[3], 2))), 2);
    }
    erv2.push_back(sqrt((sbin - 1)*sv2/sbin));
    erv3.push_back(sqrt((sbin - 1)*sv2/sbin));
    ersc32.push_back(sqrt((sbin - 1)*ssc32/sbin));
  }
  cout << "starting3" << endl;

  string name = name2 + "/obsb3_" + chg + "/ultracen_" + pname + ".dat";
  cout << name << endl;
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name << endl;
 		exit(1);
 	}


  for(int i = 0; i < olen; i++)
	{
    PRINT << M[i]/mM << " " << v2[i]/sqrt(mv2)  << " " << erv2[i] << " "<< v3[i]/sqrt(mv3) << " " << erv3[i] <<   " " <<  sc32[i]/(pow(v2[i], 2)*pow(v3[i], 2)) - 1 << " " << ersc32[i] << endl;
  }

  cout << "0-1% average multiplicity = " << mM << endl;

  PRINT.close();
}

void ucprintCGC(string name2, string chg, string pname, vector<eccs> ssort)
{
  int slen = ssort.size();

  int sbin = slen/20;

  int sc = 0;
  int npc = 0;
  vector<double> v2, v3, sc32, M, v24v22, v34v32, cn2, cn3;
  double mv2 = 0, mv3 = 0, msc32 = 0, mM = 0;

  double amean2 = 0, amean3 = 0;
  vector <double> norm, Inot;
  double I0 = 0;
  while(sc < slen)
	{
    double c2 = 0, c3 = 0, c32 = 0, cM = 0;
    double mean2=0,mean3=0;

    norm.push_back(0);
    Inot.push_back(0);
    int cur = norm.size() - 1;
    double I0i = 0;
    for(int i = sc; i < (sc + sbin); i++)
		{
      if (i >= slen) break;
      c2 += pow(ssort[i].flow_harmonics_no_decays[2], 2);
      c3 += pow(ssort[i].flow_harmonics_no_decays[3], 2);
      c32 += ssort[i].flow_harmonics_no_decays[2]*ssort[i].flow_harmonics_no_decays[3];
      mean2 += ssort[i].flow_harmonics_no_decays[2];


      mean3 += ssort[i].flow_harmonics_no_decays[3];
      cM += ssort[i].multiplicity;
      I0i += ssort[i].multiplicity*ssort[i].multiplicity;
      norm[cur]++;

    }
    mv2 += mean2;
    mv3 += mean3;
    msc32 += c32;
    mM += cM;
    amean2 += mean2;
    amean3 += mean3;
    v2.push_back(sqrt(mean2/pow(I0i, 1)));
    v3.push_back(sqrt(mean3/pow(I0i, 1)));
    sc32.push_back(c32/pow(I0i, 2));
    M.push_back(cM/norm[cur]);
    cn2.push_back(c2);
    cn3.push_back(c3);
    Inot[cur] = I0i;
    I0 += I0i;

    double s2 = -2.*(c2/pow(I0i, 2) - pow(mean2/pow(I0i, 1), 1))/pow(mean2/pow(I0i, 1), 2);
    v24v22.push_back(s2);
    double s3 = -2.*(c3/pow(I0i, 2) - pow(mean3/pow(I0i, 1), 1))/pow(mean3/pow(I0i, 1), 2);
    v34v32.push_back(s3);
    cout << v2[cur] << " " <<norm[cur] << " "<< sqrt(mean2/pow(I0i, 1)) << " " << sqrt(mean3/pow(I0i, 2)) << " " << s2 << endl;
  	sc += norm[cur];
  }




  mv2 /= pow(I0, 1);
  mv3 /= pow(I0, 1);
  amean2 /= pow(I0, 1);
  amean3 /= pow(I0, 1);
  msc32 /= pow(I0, 2)*amean2*amean3;
  mM /= slen;

	cout << "averaged" << endl;


  int olen = v2.size();
  vector<double> erv2, erv3, ersc32, errv24v22, errv34v32;
  for(int i = 0; i < olen; i++)
	{
    double sv2 = 0, sv3 = 0, ssc32 = 0;
    double vs[4] = {0}, v4s[4] = {0}, ms[4] = {0};
    for(int j = i*sbin; j < ((i + 1)*sbin); j++)
		{
      if (j >= slen) break;
      double Msub = Inot[i] - pow(ssort[j].multiplicity, 1);
      sv2 += pow(v2[i]/sqrt(mv2) - sqrt((pow(v2[i], 2)*pow(Inot[i], 1) - ssort[j].flow_harmonics_no_decays[2])/pow(Msub, 1))/sqrt(mv2), 2);
      sv3 += pow(v3[i]/sqrt(mv3) - sqrt((pow(v3[i], 2)*pow(Inot[i], 1) - ssort[j].flow_harmonics_no_decays[3])/pow(Msub, 1))/sqrt(mv3), 2);
      ssc32 += pow(sc32[i]/(pow(v2[i], 2)*pow(v3[i], 2)) - (sc32[i]*pow(Inot[i], 2) - ssort[i].flow_harmonics_no_decays[2]*ssort[i].flow_harmonics_no_decays[3])/((pow(v2[i], 2)*pow(Inot[i], 1) - ssort[j].flow_harmonics_no_decays[2])*(pow(v3[i], 2)*pow(Inot[i], 1) - ssort[j].flow_harmonic[3])), 2);



       vs[2] = cn2[i] - ssort[j].flow_harmonics_no_decays[2]*ssort[j].flow_harmonics_no_decays[2];
       ms[2] = pow(v2[i], 2)*pow(Inot[i], 1) - ssort[j].flow_harmonics_no_decays[2];
       vs[3] = cn3[i] - ssort[j].flow_harmonics_no_decays[3]*ssort[j].flow_harmonics_no_decays[3];
       ms[3] = pow(v3[i], 2)*pow(Inot[i], 1) - ssort[j].flow_harmonics_no_decays[3];

       double sub = -2.*(vs[2]/pow(Msub, 1) - pow(ms[2]/pow(Msub, 1), 2))/(pow(ms[2]/pow(Msub, 1), 2));
       v4s[2] += pow(v24v22[i] - sub, 2);
       sub = -2.*(vs[3]/pow(Msub, 1) - pow(ms[3]/pow(Msub, 2), 2))/(pow(ms[3]/pow(Msub, 1), 2));
       v4s[3] += pow(v34v32[i] - sub, 2);
    }
    erv2.push_back(sqrt((norm[i] - 1)*sv2/norm[i]));
    erv3.push_back(sqrt((norm[i] - 1)*sv2/norm[i]));
    ersc32.push_back(sqrt((norm[i] - 1)*ssc32/norm[i]));
    errv24v22.push_back(sqrt((norm[i] - 1)*v4s[2]/norm[i]));
    errv34v32.push_back(sqrt((norm[i] - 1)*v4s[3]/norm[i]));
  }

  cout << "done" << endl;

  string name = name2 + "/obsb3_" + chg + "/ultracen_" + pname + "CGC.dat";
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name << endl;
 		exit(1);
 	}


  for(int i = 0; i < olen; i++)
	{
    PRINT << M[i]/mM << " " << v2[i]/sqrt(mv2)  << " " << erv2[i] << " "<< v3[i]/sqrt(mv3) << " " << erv3[i] <<   " " <<  sc32[i]/(pow(v2[i], 2)*pow(v3[i], 2)) - 1 << " " << ersc32[i] << " "  << v24v22[i] << " " << errv24v22[i] << " "   << v34v32[i] << " "  << errv34v32[i] <<  endl;
  }

  PRINT.close();
  cout << "printed" << endl;
}


//
//
// void SortEccentricitiesIntoCentralityBins(string folder, string sort_file, vector<eccs> all,  vector< vector<eccs> > & out, vector<double> & cens,int bin){
//
// //      calculates max and min of each centrality window according to binning
// 	vector <double> clist;
// 	vector<double> cmax,cmin;
//
//
// 	for (int i=0;i<(100-bin);i+=bin) {
//
// 	cmax.push_back(i+bin);
// 	cmin.push_back(i);
// 	cens.push_back(i+bin/2.);
// 	}
//
//   int censize=cens.size();
// 	vector<double> nph,npl;
// 	nph.resize(censize);
// 	npl.resize(censize);
//
//
// 	string named="trento/"+folder+"/"+sort_file+".dat";
// 	ifstream inputd(named.c_str());
// 	if (!inputd.is_open())
//  	{
//  	cout << "Can't open " << named << endl;
//  	exit(1);
//  	}
//  	string line;
//  	int imin=0,imax=0;
//
// 	while (getline(inputd,line)) {
//
// 	vector<double> y (2,0) ;
// 	vector<double> y2 (4,0) ;
// 	vector<string> x = split(line, ' ');
// 	if (x.size()<4) vector<string> x = split(line, '\t');
//
//
//
// 	for(int j=0;j<2;j++)
// 	{
// 	stringstream s;
// 	s << x[j];
// 	s >> y[j];
// 	}
// 	for(int j=2;j<4;j++)
// 	{
// 	stringstream s;
// 	s << x[j];
// 	s >> y2[j];
// 	}
//
// 	if (cmin[imin]==y[0]){//min of the centrality bin
// 	 nph[imin]=y2[2]; //max entropy for centrality bin
// 	 imin++;
// 	 if (imin>censize) break;
// 	 }
// 	if (cmax[imax]==y[1]){ //max of the centrality bin
// 	 npl[imax]=y2[3];//min entropy for centrality bin
// 	 imax++;
// 	 if (imax>censize) break;
// 	 }
//
//
//
// 	}
// 	inputd.close();
//
//
//
// 	//sets up variables to the approparite number of bins in centrality
// 	int tot=all.size();
//
// 	out.resize(censize);
// 	for (int i=0;i<tot;i++){
// 		for (int c=0;c<censize;c++){
// 			if (all[i].si<=nph[c]&&all[i].si>npl[c]){
// 			out[c].push_back(all[i]);
// 			break;
// 			}
// 		}
// 	}
//
//
// }


double pred(double k1, double k2, double ec)
{
	return (k1*ec + k2*pow(ec, 3));
}


void CalculateCummulants(string name2, string chg, vector<vector<eccs>> & sorted_eccentricities, vector<double> & cens, int bin)
{
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//	Open Files
	//string name3="trento/"+name2+"/ND/vns1.dat";
	string name3 = name2 + "/obsb3_" + chg + "/vns1bin" + convertInt(bin) + "rat.dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name3 << endl;
 		exit(1);
 	}

  string namet4 = name2 + "/obsb3_" + chg + "/vns1bin" + convertInt(bin) + "to4.dat";
	ofstream PRINTt4(namet4.c_str());
	if (!PRINTt4.is_open())
 	{
 		cout << "Can't open " << namet4 << endl;
 		exit(1);
 	}
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	int ssize = sorted_eccentricities.size();

	for (int cent = 0; cent < ssize; cent++)
	{
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate Values
		int evs = sorted_eccentricities[cent].size();
		double avg_Vn[7] = {0}, avg_Vn_2nd[7] = {0}, avg_Vn_3rd[7] = {0};
		double Vn_4part[7] = {0}, Vn_6part[7] = {0};
		for(int ev = 0; ev < evs; ev++)
		{
			//	Calculate average Vn, Vn^2, Vn^3
			for (int n = 2; n <= 5; n++)
			{
		 		avg_Vn[7] += pow(sorted_eccentricities[cent][ev].eccentricity[n], 2)/evs;
		 		avg_Vn_2nd[n] += pow(sorted_eccentricities[cent][ev].eccentricity[n], 4)/evs;
		 		avg_Vn_3rd[n] += pow(sorted_eccentricities[cent][ev].eccentricity[n], 6)/evs;
		 	}
		}
		//	Calculate 4 and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
			Vn_4part[n] = pow(2.*pow(avg_Vn[7], 2)
												- avg_Vn_2nd[n], 0.25)
										/sqrt(avg_Vn[7]);
			Vn_6part[n] = pow(0.25*(avg_Vn_3rd[n]
													- 9.*avg_Vn[7]*avg_Vn_2nd[n]
													+ 12.*pow(avg_Vn[7], 3)), 1./6.)
										/pow(2.*pow(avg_Vn[7], 2) - avg_Vn_2nd[n], 0.25);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate Errors
		double EvErr_avg_Vn[7] = {0}, EvErr_avg_Vn_2nd[7] = {0},  EvErr_avg_Vn_3rd[7] = {0};
		double error_Vn_4part[7] = {0}, error_Vn_4part_4th[7] = {0},  error_Vn_6part[7] = {0};
		for(int ev = 0; ev < evs; ev++)
		{
			// Calculate error by event of average Vn, Vn^2, Vn^3
			for (int n = 2; n <= 5; n++)
			{
		 		EvErr_avg_Vn[n] = (avg_Vn[7]*evs - pow(sorted_eccentricities[cent][ev].eccentricity[n], 2))/(evs - 1);
		 		EvErr_avg_Vn_2nd[n] = (avg_Vn_2nd[n]*evs - pow(sorted_eccentricities[cent][ev].eccentricity[n], 4))/(evs - 1);
		 		EvErr_avg_Vn_3rd[n] = (avg_Vn_3rd[n]*evs - pow(sorted_eccentricities[cent][ev].eccentricity[n], 6))/(evs - 1);
		 	}

			//	Calculate error for 4 and 6 particle cummulants (Also 4 particle cummulant to the 4th power)
		 	for (int n = 2; n <= 5; n++)
			{
		 		error_Vn_4part[n] += pow(Vn_4part[n]
																- pow(2.*pow(EvErr_avg_Vn[n], 2)
																 			- EvErr_avg_Vn_2nd[n], 0.25)
														/sqrt(EvErr_avg_Vn[n]), 2);
      	error_Vn_4part_4th[n] += pow(pow(Vn_4part[n], 4)
																		- (2.*pow(EvErr_avg_Vn[n], 2)
																				- EvErr_avg_Vn_2nd[n])
																	/pow(EvErr_avg_Vn[n], 2), 2);
				error_Vn_6part[n] += pow(Vn_6part[n]
																- pow(0.25*(EvErr_avg_Vn_3rd[n]
																			- 9.*EvErr_avg_Vn[n]*EvErr_avg_Vn_2nd[n]
																			+ 12.*pow(EvErr_avg_Vn[n], 3)), 1./6.)
														/pow(2.*pow(EvErr_avg_Vn[n], 2)
																- EvErr_avg_Vn_2nd[n], 0.25), 2);
			}
		}
		//	Normalize Errors
		for (int n = 2; n <= 5; n++)
		{
			error_Vn_4part[n] = sqrt(error_Vn_4part[n]*(evs - 1)/evs);
      error_Vn_4part_4th[n] = sqrt(error_Vn_4part_4th[n]*(evs - 1)/evs);
			error_Vn_6part[n] = sqrt(error_Vn_6part[n]*(evs - 1)/evs);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Output calculated values and errors
		PRINT << cens[cent] << " " ;

		for (int n = 2; n <= 5; n++) PRINT <<  Vn_4part[n] << " " << error_Vn_4part[n] << " ";
		for (int n = 2;n <= 5; n++) PRINT <<  Vn_6part[n] << " " << error_Vn_6part[n] << " ";
		PRINT << endl;

    PRINTt4 << cens[cent] << " " ;

    for (int n = 2;n <= 5; n++) PRINTt4 <<  pow(Vn_4part[n], 4) << " " << error_Vn_4part_4th[n] << " ";
    PRINTt4 << endl;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	}

	PRINT.close();
  PRINTt4.close();
}


void radii(string name2, string chg, vector<vector<eccs>> & sorted_eccentricities, vector<double> & cens, int bin)
{


	string name3 = name2 + "/obsb3_" + chg + "/radii_" + convertInt(bin) + ".dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name3 << endl;
 		exit(1);
 	}

	int ssize = sorted_eccentricities.size();

	for (int j = 0; j < ssize; j++)
	{
		int evs = sorted_eccentricities[j].size();
		double r2 = 0;
		for(int ev = 0; ev < evs; ev++)
		{
		  r2 += sorted_eccentricities[j][ev].radius*sorted_eccentricities[j][ev].radius;
		}
		r2 /= evs;
		PRINT << cens[j] << " " << sqrt(r2) << endl;
	}

	PRINT.close();
}

void vCGC(string name2, string chg, vector<vector<eccs>> & sorted_eccentricities, vector<double> & cens, int bin)
{

	//string name3="trento/"+name2+"/ND/vns1.dat";
	string name3 = name2 + "/obsb3_" + chg + "/CGCfluc" + convertInt(bin) + ".dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name3 << endl;
 		exit(1);
 	}

	int ssize = sorted_eccentricities.size();

	double v2cen = 0, v3cen = 0, Icen = 0;
	for (int j = 0; j < ssize; j++)
	{
		int evs = sorted_eccentricities[j].size();
		double v[4] = {0}, v4[4] = {0}, vnm = 0, err[4] = {0}, mv[4] = {0}, v2[4] = {0};
		double I0i = 0;
		for(int ev = 0; ev < evs; ev++)
		{
			I0i += sorted_eccentricities[j][ev].multiplicity*sorted_eccentricities[j][ev].multiplicity;
		 	for (int n = 2; n <= 3; n++)
			{
      	mv[n] += sorted_eccentricities[j][ev].flow_harmonics_no_decays[n];
		 		v[n] += sorted_eccentricities[j][ev].flow_harmonics_no_decays[n]*sorted_eccentricities[j][ev].flow_harmonics_no_decays[n];
		 	}
    	vnm += sorted_eccentricities[j][ev].flow_harmonics_no_decays[2]*sorted_eccentricities[j][ev].flow_harmonics_no_decays[3];
		}
		for (int n = 2; n <= 3; n++) v4[n] = -2.*(v[n]/evs - pow(mv[n]/evs, 2))/pow(mv[n]/evs, 2);

    vnm /= evs;
    double sc32 = vnm/(mv[2]/evs*mv[3]/evs) - 1;
    for (int n = 2; n <= 3; n++) v2[n] = sqrt(mv[n]/I0i);
    if (j == 0)
		{
    	v2cen = v2[2];
    	v3cen = v2[3];
    	Icen = I0i;
    }


		double vs[4] = {0}, v4s[4] = {0}, ms[4] = {0};
    double ssc32 = 0;
		for(int ev = 0; ev < evs; ev++)
		{
			for (int n = 2; n <= 3; n++)
			{
		 		vs[n] = v[n] - sorted_eccentricities[j][ev].flow_harmonics_no_decays[n]*sorted_eccentricities[j][ev].flow_harmonics_no_decays[n];
        ms[n] = mv[n] - sorted_eccentricities[j][ev].flow_harmonics_no_decays[n];
		 	}

		 	for (int n = 2; n <= 3; n++)
			{
       	double sub = -2.*(vs[n]/(evs - 1) - pow(ms[n]/(evs - 1), 2))/(pow(ms[n]/(evs - 1), 2));
		 		v4s[n] += pow(v4[n] - sub, 2);
			}
      ssc32 += pow(sc32 + 1 - (evs - 1)*(vnm*evs - sorted_eccentricities[j][ev].flow_harmonics_no_decays[2]*sorted_eccentricities[j][ev].flow_harmonics_no_decays[3])/((vs[2]*evs - sorted_eccentricities[j][ev].flow_harmonics_no_decays[2])*(vs[3]*evs - sorted_eccentricities[j][ev].flow_harmonics_no_decays[3])), 2);

		}
		for (int n = 2; n <= 3; n++)
		{
			err[n] = sqrt((evs - 1)*v4s[n]/evs);
		}
    err[0] = sqrt((evs - 1)*ssc32/evs);
		PRINT << cens[j] << " " << I0i/Icen << " " <<  v2[2]/v2cen <<  " "  <<  v2[3]/v3cen <<  " ";
		for (int n = 2; n <= 3; n++) PRINT <<  v4[n] << " " << err[n] << " ";
		PRINT << sc32 << " " << err[0] <<  endl;


	}

	PRINT.close();
}
//void vnsb1(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens){

//	//string name3="trento/"+name2+"/ND/vns1.dat";
//	string name3="trento/"+name2+"/obsb3_"+chg+"/vns1bin1n.dat";
//	ofstream PRINT(name3.c_str());
//	if (!PRINT.is_open())
// 	{
// 	cout << "Can't open " << name3 << endl;
// 	exit(1);
// 	}

//	int ssize=sorted_eccentricities.size();
//
//	for (int j=0;j<ssize;j++ ){
//		int evs=sorted_eccentricities[j].size();
//		double v[7]={0},v4[7]={0},v6[7]={0},sv4[7]={0},sv6[7]={0},err[7]={0},err6[7]={0};
//		for(int ev=0;ev<evs;ev++){
//		 for (int n=2;n<=5;n++) {
//		 	v[n]+=sorted_eccentricities[j][ev].eccentricity[n]*sorted_eccentricities[j][ev].eccentricity[n];
//		 	sv4[n]+=pow(sorted_eccentricities[j][ev].eccentricity[n],4);
//		 	sv6[n]+=pow(sorted_eccentricities[j][ev].eccentricity[n],6);
//		 }
//		}
//		for (int n=2;n<=5;n++) {
//		v4[n]=pow(2.*v[n]/evs*v[n]/evs-sv4[n]/evs,0.25)/sqrt(v[n]/evs);
//		v6[n]=pow(0.25*(sv6[n]/evs-9.*v[n]/evs*sv4[n]/evs+12.*pow(v[n]/evs,3) ),1./6.)/pow(2.*v[n]/evs*v[n]/evs-sv4[n]/evs,0.25);
//		}
//
//
//		double vs[7]={0},v4s[7]={0},sv4s[7]={0},v6s[7]={0},sv6s[7]={0};
//		for(int ev=0;ev<evs;ev++){
//		 for (int n=2;n<=5;n++) {
//		 	vs[n]=v[n]-sorted_eccentricities[j][ev].eccentricity[n]*sorted_eccentricities[j][ev].eccentricity[n];
//		 	sv4s[n]=sv4[n]-pow(sorted_eccentricities[j][ev].eccentricity[n],4);
//		 	sv6s[n]=sv6[n]-pow(sorted_eccentricities[j][ev].eccentricity[n],6);
//		 }
//
//		 for (int n=2;n<=5;n++) {
//		 	v4s[n]+=pow(v4[n]-pow(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1),0.25)/sqrt(vs[n]/(evs-1)),2);
//			v6s[n]+=pow(v6[n]-pow(0.25*(sv6s[n]/(evs-1)-9.*vs[n]/(evs-1)*sv4s[n]/(evs-1)+12.*pow(vs[n]/(evs-1),3) ),1./6.)/pow(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1),0.25),2);
//
//
//			}
//
//		}
//		for (int n=2;n<=5;n++) {
//			err[n]=sqrt((evs-1)*v4s[n]/evs);
//			err6[n]=sqrt((evs-1)*v6s[n]/evs);
//			}
//
//		PRINT << cens[j] << " " ;
//
//		for (int n=2;n<=5;n++) PRINT <<  v4[n] << " " << err[n] << " ";
//		for (int n=2;n<=5;n++) PRINT <<  v6[n] << " " << err6[n] << " ";
//		PRINT << endl;
//
//
//	}
//
//	PRINT.close();
//}



void CalculateCummulantsWithoutRatios(string name2, string chg, vector<vector<eccs>> & sorted_eccentricities, vector<double> & cens, int bin)
{

//	string name3="trento/"+name2+"/ND/vns1bin5.dat";
//	string name4="trento/"+name2+"/ND/vn2err.dat";
	string name3 = name2 + "/obsb3_" + chg + "/enstotbin" + convertInt(bin) + ".dat";
	string name4 = name2 + "/obsb3_" + chg + "/en2errbin" + convertInt(bin) + ".dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 		cout << "Can't open " << name3 << endl;
 		exit(1);
 	}



	ofstream PRINT2(name4.c_str());
	if (!PRINT2.is_open())
 	{
 		cout << "Can't open " << name4 << endl;
 		exit(1);
 	}

	int ssize = sorted_eccentricities.size();

	for (int cent = 0; cent < ssize; cent++)
	{
		int evs = sorted_eccentricities[cent].size();
		double avg_Vn[7] = {0}, avg_Vn_2nd[7] = {0}, avg_Vn_3rd[7] = {0};
		double Vn_2part[7] = {0}, Vn_4part[7] = {0}, Vn_6part[7] = {0};
		double V2_2part_V3_2part = 0;

		double mean_pt = 0, mean_pt_2nd = 0, RMS_pt = 0;

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values
		for(int ev = 0; ev < evs; ev++) mean_pt += sorted_eccentricities[cent][ev].average_pt;
		mean_pt /= evs;

		for(int ev = 0; ev < evs; ev++)
		{
			//	Calculate average Vn, Vn^2, Vn^3
			for (int n = 2; n <= 5; n++)
			{
		 		avg_Vn[n] += pow(sorted_eccentricities[cent][ev].eccentricity[n], 2)/evs;
		 		avg_Vn_2nd[n] += pow(sorted_eccentricities[cent][ev].eccentricity[n], 4)/evs;
		 		avg_Vn_3rd[n] += pow(sorted_eccentricities[cent][ev].eccentricity[n], 6)/evs;
		 	}

		 	mean_pt_2nd += pow(sorted_eccentricities[cent][ev].average_pt - mean_pt, 2.);
		}
		//	Calculate 2, 4, and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
			Vn_2part[n] = sqrt(avg_Vn[n]);
			Vn_4part[n] = 2.*pow(avg_Vn[n], 2) - avg_Vn_2nd[n];
			Vn_6part[n] = 0.25*(avg_Vn_3rd[n] - 9.*avg_Vn[n]*avg_Vn_2nd[n] + 12.*pow(avg_Vn[n], 3));
		}
		V2_2part_V3_2part = sqrt(avg_Vn[2]/avg_Vn[3]);

		mean_pt_2nd /= evs;
		RMS_pt = sqrt(mean_pt_2nd)/mean_pt;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
		double EvErr_avg_Vn[7] = {0}, EvErr_avg_Vn_2nd[7] = {0},  EvErr_avg_Vn_3rd[7] = {0};
		double err_Vn_2part[7] = {0}, err_Vn_4part[7] = {0},  err_Vn_6part[7] = {0};
		double err_V2_2part_V3_2part = 0;

		double err_mean_pt = 0, err_mean_pt_2nd = 0, err_RMS_pt = 0;

		for(int ev = 0; ev < evs; ev++)
		{
			// Calculate error by event of average Vn, Vn^2, Vn^3
			for (int n = 1; n <= 6; n++)
			{
		 		EvErr_avg_Vn[n] = (avg_Vn[n]*evs - pow(sorted_eccentricities[cent][ev].eccentricity[n], 2))/(evs - 1);
		 		EvErr_avg_Vn_2nd[n] = (avg_Vn_2nd[n]*evs - pow(sorted_eccentricities[cent][ev].eccentricity[n], 4))/(evs - 1);
		 		EvErr_avg_Vn_3rd[n] = (avg_Vn_3rd[n]*evs - pow(sorted_eccentricities[cent][ev].eccentricity[n], 6))/(evs - 1);
		 	}

			//	Calculate error for 2, 4, and 6 particle cummulants
		 	for (int n = 2; n <= 5; n++)
			{
		 		err_Vn_2part[n] += pow(Vn_2part[n] - sqrt(EvErr_avg_Vn[n]), 2);
				err_Vn_4part[n] += pow(Vn_4part[n] - (2.*pow(EvErr_avg_Vn[n], 2) - EvErr_avg_Vn_2nd[n]), 2);
				err_Vn_6part[n] += pow(Vn_6part[n] - 0.25*(EvErr_avg_Vn_3rd[n] - 9.*EvErr_avg_Vn[n]*EvErr_avg_Vn_2nd[n] + 12.*pow(EvErr_avg_Vn[n], 3)), 2);
			}

			err_V2_2part_V3_2part += pow(V2_2part_V3_2part - sqrt(EvErr_avg_Vn[2]/EvErr_avg_Vn[3]), 2);

			double EvErr_mean_pt = (mean_pt*evs - sorted_eccentricities[cent][ev].average_pt)/(evs - 1.);
			double EvErr_mean_pt_2nd = (mean_pt_2nd*evs - pow(sorted_eccentricities[cent][ev].average_pt - subpt, 2))/(evs - 1);

		 	err_mean_pt += pow(mean_pt - EvErr_mean_pt, 2);
		 	err_mean_pt_2nd += pow(mean_pt_2nd - EvErr_mean_pt_2nd, 2);
		 	err_RMS_pt = pow(RMS_pt - sqrt(EvErr_mean_pt_2nd)/EvErr_mean_pt, 2);
		}

		//	Normalize errors
		for (int n = 2; n <= 5; n++)
		{
			err_Vn_2part[n] = sqrt(err_Vn_2part[n]*(evs - 1)/evs);
			err_Vn_4part[n] = sqrt(err_Vn_4part[n]*(evs - 1)/evs);
			err_Vn_6part[n] = sqrt(err_Vn_6part[n]*(evs - 1)/evs);
		}
		err_V2_2part_V3_2part = sqrt(err_V2_2part_V3_2part*(evs - 1)/evs);

		err_mean_pt = sqrt(err_mean_pt*(evs - 1)/evs);
		err_mean_pt_2nd = sqrt(err_mean_pt_2nd*(evs - 1)/evs);
		err_RMS_pt = sqrt(err_RMS_pt*(evs - 1)/evs);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Output calculated values and errors
		PRINT << cens[cent] << " " ;

		for (int n = 2; n <= 5; n++) PRINT <<  Vn_4part[n] << " " << err_Vn_4part[n] << " ";
		for (int n = 2; n <= 5; n++) PRINT <<  Vn_6part[n] << " " << err_Vn_6part[n] << " ";
		PRINT << V2_2part_V3_2part << " " << err_V2_2part_V3_2part << endl;

		PRINT2 << cens[cent] << " " ;
		for (int n = 2; n <= 5; n++) PRINT2 <<  Vn_2part[n] << " " << err_Vn_2part[n] << " ";
		PRINT2 << endl;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	}

	PRINT.close();
	PRINT2.close();
}


void PredictCummulants(string full, string runtype,  string runname,string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens, int bin)
{


	// Reads in the previously calculated kappa's for the linear+cupic response
  string name = runtype + "/obschgv1/nonlin_bin" + convertInt(bin) + ".dat";
  ifstream input(name.c_str());
  if (!input.is_open())
  {
    cout << "Can't open predictions file " << name << endl;
    return;
  }

  string line;
  int cc = 0;
  vector <double> k1v2, k2v2;
  vector <double> k1v3, k2v3;
  while (getline(input, line))
	{
    vector<double> y (13, 0) ;
    vector<string> x = split(line, ' ');


    for(int j = 0; j < 13; j++)
    {
    	stringstream s;
    	s << x[j];
    	s >> y[j];
    }
    k1v2.push_back(y[4 - 1]);
    k2v2.push_back(y[6 - 1]);
    k1v3.push_back(y[10 - 1]);
    k2v3.push_back(y[12 - 1]);
    cc++;
	}
  input.close();
	// Done reading in kappas



	// define names of outputfiles
  string name6 = full + "/obsb3_" + chg + "/predictv24v22_bin" + convertInt(bin) + "_" + runname + ".dat";

  ofstream PRINT4(name6.c_str());
  if (!PRINT4.is_open())
  {
  	cout << "Can't open " << name6 << endl;
  	exit(1);
  }

  string name2 = full + "/obsb3_" + chg + "/predictv22_bin" + convertInt(bin) + "_" + runname + ".dat";

  ofstream PRINT(name2.c_str());
  if (!PRINT.is_open())
  {
  	cout << "Can't open " << name2 << endl;
  	exit(1);
  }
  //end creating outputfiles

//calculates vn{2} and vn{4}/vn{2} using linear+cubic response from the eccentricities (no actual flow used here)
// j is the centrality bin
	int ssize = sorted_eccentricities.size();
	for (int j = 0; j < ssize; j++)
	{
		int evs = sorted_eccentricities[j].size();
		double err[7] = {0}, err2[7] = {0};

		//predicts v2{4}/v2{2} from eccentricities
		double ecc22[4] = {0}, ecc44[4] = {0}, pecv2[4] = {0}, pecv24[4] = {0};
		for(int ev = 0; ev < evs; ev++)
		{
		 	for (int n = 2; n <= 3; n++)
		 	{
      	double en = 0;
		 		if (n == 2)en = pred(k1v2[j],k2v2[j], sorted_eccentricities[j][ev].eccentricity[n]);
      	if (n == 3)en = pred(k1v3[j],k2v3[j], sorted_eccentricities[j][ev].eccentricity[n]);
		 		ecc22[n] += pow(en, 2);
		 		ecc44[n] += pow(en, 4);
		 	}// calculates moments of the distributions
		}
		for (int n = 2; n <= 3; n++)
		{
			ecc22[n] /= evs;
			ecc44[n] /= evs;
      pecv2[n] = sqrt(ecc22[n]);
			pecv24[n] = pow(2*pow(ecc22[n], 2) - ecc44[n], 0.25)/sqrt(ecc22[n]);
		} // these are the central values, error bars calculated below


   //here we used jackknife error to estimate the error bands
		double vs[7] = {0}, v2s[7] = {0}, v4s[7] = {0}, sv4s[7] = {0};
		for(int ev = 0; ev < evs; ev++)
		{
		 	for (int n = 2; n <= 3; n++)
			{
      	double en = 0;
      	if (n == 2)	en = pred(k1v2[j],k2v2[j], sorted_eccentricities[j][ev].eccentricity[n]);
      	if (n == 3)	en = pred(k1v3[j],k2v3[j], sorted_eccentricities[j][ev].eccentricity[n]);
		 		vs[n] = (ecc22[n]*evs - pow(en, 2))/(evs - 1);
		 		sv4s[n] = (ecc22[n]*evs - pow(en, 4))/(evs - 1);
		 	}

		 	for (int n = 2; n <= 3; n++)
			{
		 		v2s[n] += pow(pecv2[n] - sqrt(vs[n]), 2);
				v4s[n] += pow(pecv24[n] - pow(2.*vs[n]*vs[n] - sv4s[n], 0.25), 2);
    	}
		}


		for (int n = 2; n <= 3; n++)
		{
			err[n] = sqrt((evs - 1)*v2s[n]/evs);
			err2[n] = sqrt((evs - 1)*v4s[n]/evs);
		}
    //done calculating errors


    //printing off vn{4} and vm{2} in separate files


		PRINT4 << cens[j] << " " ;
		for (int n = 2; n <= 3; n++) PRINT4 << pecv24[n]  << " " << err2[n] << " " ;
		PRINT4 << endl;


    PRINT << cens[j] << " " ;
    for (int n = 2; n <= 3; n++) PRINT << pecv2[n]  << " " << err[n] << " " ;
    PRINT << endl;
  }


  PRINT.close();
  PRINT4.close();


}


//void vns3(string name2, string chg, vector< vector<eccs> > & sorted_eccentricities, vector<double> & cens){

////	string name3="trento/"+name2+"/ND/vns1bin5.dat";
////	string name4="trento/"+name2+"/ND/vn2err.dat";
//	string name3="trento/"+name2+"/obsb3_"+chg+"/vns1bin5.dat";
//	string name4="trento/"+name2+"/obsb3_"+chg+"/vn2err.dat";
//	ofstream PRINT(name3.c_str());
//	if (!PRINT.is_open())
// 	{
// 	cout << "Can't open " << name3 << endl;
// 	exit(1);
// 	}
//
//
//
//	ofstream PRINT2(name4.c_str());
//	if (!PRINT2.is_open())
// 	{
// 	cout << "Can't open " << name4 << endl;
// 	exit(1);
// 	}

//	int ssize=sorted_eccentricities.size();
//
//	for (int j=0;j<ssize;j++ ){
//		int evs=sorted_eccentricities[j].size();
//		double v[7]={0},v2[7]={0},v4[7]={0},v6[7]={0},sv4[7]={0},sv6[7]={0},err2[7]={0},err[7]={0},err6[7]={0};
//		double mpt=0, mpt2=0,rpt=0;
//		double v2v3=0;
//
//		for(int ev=0;ev<evs;ev++) mpt+=sorted_eccentricities[j][ev].average_pt;
//		mpt/=evs;
//
//		for(int ev=0;ev<evs;ev++){
//		 for (int n=2;n<=5;n++) {
//		 	v[n]+=pow(sorted_eccentricities[j][ev].eccentricity[n],2);
//		 	sv4[n]+=pow(sorted_eccentricities[j][ev].eccentricity[n],4);
//		 	sv6[n]+=pow(sorted_eccentricities[j][ev].eccentricity[n],6);
//		 }
//
//		 mpt2+=pow(sorted_eccentricities[j][ev].average_pt-mpt,2.);
//		}
//		for (int n=2;n<=5;n++) {
//		v2[n]=sqrt(v[n]/evs);
//		v4[n]=2.*v[n]/evs*v[n]/evs-sv4[n]/evs;
//		v6[n]=0.25*(sv6[n]/evs-9.*v[n]/evs*sv4[n]/evs+12.*pow(v[n]/evs,3) );
//		}
//		v2v3=sqrt(v[2]/v[3]);
//
//		mpt2/=evs;
//		rpt=sqrt(mpt2)/mpt;
//
//		double vs[7]={0},v2s[7]={0},v4s[7]={0},sv4s[7]={0},v6s[7]={0},sv6s[7]={0},sv2v3=0,ev2v3=0;
//		double empt=0, empt2=0,erpt=0;
//		double smpt=0, smpt2=0,srpt=0;
//		for(int ev=0;ev<evs;ev++){
//		 for (int n=1;n<=6;n++) {
//		 	vs[n]=v[n]-pow(sorted_eccentricities[j][ev].eccentricity[n],2);
//		 	sv4s[n]=sv4[n]-pow(sorted_eccentricities[j][ev].eccentricity[n],4);
//		 	sv6s[n]=sv6[n]-pow(sorted_eccentricities[j][ev].eccentricity[n],6);
//		 }
//
//		 for (int n=2;n<=5;n++) {
//		 	v2s[n]+=pow(v2[n]-sqrt(vs[n]/(evs-1)),2);
//			v4s[n]+=pow(v4[n]-(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1) ),2);
//			v6s[n]+=pow(v6[n]-0.25*(sv6s[n]/(evs-1)-9.*vs[n]/(evs-1)*sv4s[n]/(evs-1)+12.*pow(vs[n]/(evs-1),3) ),2);
//			}
//			double subpt=(mpt*evs-sorted_eccentricities[j][ev].average_pt)/(evs-1.);
//
//			double subpt2=(mpt2*evs-pow(sorted_eccentricities[j][ev].average_pt-subpt,2) )/(evs-1);
//
//			sv2v3+=pow(v2v3-sqrt(vs[2]/vs[3]),2);
//		 	smpt+=pow(mpt-subpt,2);
//		 	smpt2+=pow(mpt2-subpt2,2);
//		 	srpt=pow(rpt-sqrt(subpt2)/subpt,2);
//		}
//
//
//
//		for (int n=2;n<=5;n++) {
//			err2[n]=sqrt((evs-1)*v2s[n]/evs);
//			err[n]=sqrt((evs-1)*v4s[n]/evs);
//			err6[n]=sqrt((evs-1)*v6s[n]/evs);
//		}
//		ev2v3=sqrt((evs-1)*sv2v3/evs);
//		empt=sqrt((evs-1)*smpt/evs);
//		empt2=sqrt((evs-1)*smpt2/evs);
//		erpt=sqrt((evs-1)*srpt/evs);
//
//		PRINT << cens[j] << " " ;
//
//		for (int n=2;n<=5;n++) PRINT <<  v4[n] << " " << err[n] << " ";
//		for (int n=2;n<=5;n++) PRINT <<  v6[n] << " " << err6[n] << " ";
//		PRINT << v2v3 << " " << ev2v3 << endl;
//
//		PRINT2 << cens[j] << " " ;
//		for (int n=2;n<=5;n++) PRINT2 <<  v2[n] << " " << err2[n] << " ";
//		PRINT2 << endl;
//	}
//
//	PRINT.close();
//	PRINT2.close();
//}

//*******************************************************************************************************************************************************************************
//End of function definitions
//*******************************************************************************************************************************************************************************
