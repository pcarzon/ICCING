#include <iostream>
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

// sort - sorting events into centrality bins (actually does the work) (sall[centrality][i])
void sort(string name,string sfile,std::vector<eccs> all,  std::vector< std::vector<eccs> > & out, std::vector<double> & cens,int bin);

// NPsc - symetric cummulants without fancy stuff as one does in mathematica
void NPsc(string name2,  string chg,pars p,  std::vector<eccs>  & ec);

// bsort - 
void bsort(string folder, std::vector<eccs> all);

// maxnp - 
int maxnp( std::vector<eccs> s1)
{
	int ss=s1.size();
	int max=0;
	for (int i=0;i<ss;i++){
		if (s1[i].npart>max) max=s1[i].npart;
	}
    return max;
}

// pred - 
double pred(double k1, double k2, double ec);

// predcum - predict v22 and v24 using kappas (skanda's work) comment out for now
void predcum(string ictype,string runtype, string runname,string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens, int bin);

// sortmake - sort in centrality (make file of cuttoffs)
double sortmake(string folder, string sort_file, std::vector<eccs> all);

// sortmake2 - sort in npart (make file of cuttoffs)
int sortmake2(string folder, std::vector<eccs> all);

//void sort2(string name,std::vector<eccs> all,  std::vector< std::vector<eccs> > & out, std::vector<double> & cens);

//void sort3(string name,std::vector<eccs> all,  std::vector< std::vector<eccs> > & out, std::vector<double> & cens);

// ultracen - 
void ultracen(string name2, string chg, std::vector<eccs> all, double scut, int npcut);

// ucprintCGC - 
void ucprintCGC(string name2,string chg,string pname,vector<eccs> ssort);

// ucprint - 
void ucprint(string name2,string chg,string pname,vector<eccs> ssort);

// vCGC - 
void vCGC(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin);

// readecc - read in eccentricities from folder (add an if to skip nans)
int readecc(string name2, vector<eccs> & all);

// vns - calculate all Vnm's with error bars (easy mode with statistics)
void vns(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin);

// vnsnorat - same as vns but with out stuff like v24/v22
void vnsnorat(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin);
//void vnsb1(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens);

//void vns2(string name2,string chg,  std::vector< std::vector<eccs> > & sall, std::vector<double> & cens);

//void vns3(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens);

// radii - 
void radii(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin);


// compareByLength - compare by multiplicity
bool compareByLength(const eccs &a, const eccs &b)
{
    return a.si > b.si;
}

// compareByNpart - 
bool compareByNpart(const eccs &a, const eccs &b)
{
    return a.npart > b.npart;
}
//*******************************************************************************************************************************************************************************
//End of function declarations
//*******************************************************************************************************************************************************************************

//*******************************************************************************************************************************************************************************
//main
//*******************************************************************************************************************************************************************************
int main (int argc, char *argv[]) //execute with ./a.out folder list_name sortfile_name PID
{
	cout << "HELLO!!!!!!" << endl;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//assing input to strings for output of data  
	
	// name = folder path to file of eccentricities (assumption /trento/name)
	// sname = name of list of eccentricities without .dat extension
	// chg = set to 0
	// sfile = output file name
	// runtype = folder for running with kappas (from Skanda)
	// runname = kappa file
	string name,sname,chg,sfile,runtype,runname;
	
	if (argv[1])
	{
		name=argv[1];

	}
	else name="";

	if (argv[2])
	{
		sname=argv[2];

	}
	else sname="";

	if (argv[3])
	{
		sfile=argv[3];

	}
	else sfile="";

	if (argv[4])
	{
		chg=argv[4];

	}
	else chg="0";
	  
    if (argv[5])//specify IC folder to obtain kappas
	{
		runtype=argv[5];

	}
	else runtype="0";

    if (argv[6])//specify IC folder to obtain kappas
	{
		runname=argv[6];

	}
	else runname="0";
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	vector<eccs> all,aecc;
	int tot=0;
	string full=name+"/"+sname;
  cout << full << endl;
	int maxl=readecc(full, all);
	cout << " out maxl=" << maxl << endl;
  double scut=sortmake(name, sfile, all);
  int npcut=sortmake2(name, all);

	cout << all.size() << endl;
	vector< vector<eccs> > sall,sall2,sall3;



	vector<double> cens,cens2,cens3;
  bsort(name, all);
	sort(name,sfile,all,sall,cens,5); //5% bins


//
	int ssize=sall.size();
//	for (int j=0;j<ssize;j++ ) SC.correct(sall[j]);

	string bin="5";
	vns(name,chg,sall,cens,5);
	if (maxl>13) radii(name,chg,sall,cens,5);
	vnsnorat(name,chg, sall, cens,5);
predcum(name,runtype,runname,chg,sall,cens,5);

	cout << "sort into 1 bins" << endl;
	sort(name,sfile,all,sall2,cens2,1);
	cout << "done 1%" << endl;
cout << "will do maxl " << maxl <<  endl;
  if (maxl>8) vCGC(name,chg,sall2,cens2,1);

  vns(name,chg,sall2,cens2,1);
  if (maxl>13) radii(name,chg,sall2,cens2,1);
  vnsnorat(name,chg, sall2, cens2,1);

  if (maxl>8) ultracen(name, chg,all, scut,npcut);
//  predcum(name,runtype,runname,chg,sall2,cens2,1);

	sc SC;
	string pre=name+"/obsb3_"+chg+"/";
	string into;

	pars p1;
	p1.n=3;
	p1.m=2;
	cout << "np star" << endl;
        NPsc(name,chg,p1,all);
	into=pre+"SCM32comb";
	cout << "start sc" << endl;
	SC.runec(into,p1,sall2);
	cout << "sc 32" << endl;

//	vns3(name,chg,sall,cens);

//	string bin="1";
//	sort2(name,all,sall2,cens2); //1% bins
//
//	sc SC2;
//
//	ssize=sall2.size();
//	//for (int j=0;j<ssize;j++ ) SC2.correct(sall2[j]);

//	vnsb1(name,chg,sall2,cens2);
//	vns2(name,chg, sall2,cens2);
//
//
//	sc SC;


//	string n="obs/"+name;
//	SC.eccsout(n,sall);
//	exit(0);

//	string pre="trento/"+name+"/obsb3_"+chg+"/";
//	string into;
//
//	pars p1;
//	p1.n=3;
//	p1.m=2;
//	into=pre+"SCM32comb";
//	SC.scnobin(into,p1,sall2);
//
	p1.n=4;
	p1.m=2;
	into=pre+"SCM42comb";
	SC.runec(into,p1,sall2 );

	p1.n=4;
	p1.m=3;
	into=pre+"SCM43comb";
	SC.runec(into,p1,sall2);


}
//*******************************************************************************************************************************************************************************
//End of main
//*******************************************************************************************************************************************************************************


//*******************************************************************************************************************************************************************************
//Function Definitions
//*******************************************************************************************************************************************************************************
int readecc(string name2, vector<eccs> & all){
	cout << "Now!!!!!!" << endl;
	int maxl;
	string name=name2+".dat";
	ifstream input(name.c_str());
	if (!input.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}
//	cout << "1" << endl;
	string line;
	int cc=0;
	while (getline(input,line)) {

	if (cc>3000000) break;
//	if (line.size() < 139) break;
	std::vector<double> y (15,0) ;
	//cout << line.size() << endl;
//	line.replace(line.begin(), line.end(), "\t", " ");
	std::vector<std::string> x = split(line, ' ');
//	if (x.size() < 8) break;
//	cout << "xsize " << x.size() << endl;
        int tsize=x.size();
  	if (tsize<8) break;

    maxl=8;
  //  cout << "2" << endl;
   // if (x[10] == "-nan") maxl=8; 
	
    eccs sub;
	for(int j=0;j<maxl;j++)
	{
	
//	y[j] = stod(x[j]);
	stringstream s;
	s << x[j];
	s >> y[j];
//	cout << y[j] << " " << endl;
	}

	if (maxl==5){// old mckln files
		sub.si=y[0]; 
		sub.M=y[0];
	  sub.npart=y[0];
		sub.ec[2]=y[1];
		sub.ec[3]=y[3];
		sub.ec[4]=0;
		sub.ec[5]=0;
		sub.v[2]=y[1];
		sub.v[3]=y[3];
		sub.v[4]=0;
    sub.b=0;
	}

//	cout << "3" << endl;
	//work here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	if (maxl>5){//original trento output
 sub.b=y[1];
	sub.si=y[3];
	sub.M=y[3];
 sub.npart=y[2];
	sub.ec[2]=y[4];
	sub.ec[3]=y[5];
	sub.ec[4]=y[6];
	sub.ec[5]=y[7];
	sub.v[2]=y[4];
	sub.v[3]=y[5];
	sub.v[4]=y[6];
}

	//this is all of trento stuff
  if (maxl>8){// trento+CGCoutput
	sub.vND[1]=y[8];//equivalent to "multiplicity"
  sub.vND[2]=pow(y[8],1)*y[9];
	sub.vND[3]=pow(y[8],1)*y[10];
	sub.vND[4]=pow(y[8],1)*y[11];
	 if (maxl>13) sub.R=y[13];
	 
}
//cout << "try" << endl;
	all.push_back(sub);
	cc++;
	}
cout << "yeah" << endl;
	input.close();

  cout <<"done" << endl;

	return maxl;
}

void bsort(string folder, std::vector<eccs> all){

  double becc2[12]={0},becc3[12]={0},cnt[12]={0};

  int tot=all.size();
	for (int i=0;i<tot;i++){
		for (int b=0;b<12;b++){
      int bmin=b-1;
      if (bmin<0) bmin=0;
			if (all[i].b<=b&&all[i].b>=bmin) {
        becc2[b]+=all[i].ec[2];
        becc3[b]+=all[i].ec[3];
        cnt[b]++;
        break;
      }

		}
	}



  for (int b=0;b<12;b++){
    cout << b << " " << cnt[b] << " " << becc2[b]/cnt[b] << " " << becc3[b]/cnt[b] << endl;
  }
}


double sortmake(string folder, string sort_file, std::vector<eccs> all){

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax,cmin;


	cout << "start sort" << endl;
  std::sort (all.begin(), all.end(), compareByLength);

  int dbin= int (all.size()/100);
  //cout << dbin << " " << all.size() << endl;

	vector<double> cuts;
	for (int c=0;c<100;c++) {
    cuts.push_back(all[c*dbin].si);
  }

	int end=dbin*100-1;

	cuts[100]=all[end].si;


  string name=folder+"/"+sort_file+".dat";
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}

  for (int c=0;c<99;c++) PRINT << c << " " << c+1 << " " <<  cuts[c]  << " "  << cuts[c+1] << endl;

  PRINT.close();
  double cuts1=cuts[1];


  string name2=folder+"/"+sort_file+"5.dat";
  ofstream PRINT2(name2.c_str());
	if (!PRINT2.is_open())
 	{
 	cout << "Can't open " << name2 << endl;
 	exit(1);
 	}

  for (int c=0;c<99;c=c+5) PRINT2 << c << " " << c+5 << " " <<  cuts[c]  << " "  << cuts[c+5] << endl;

  PRINT2.close();

  double cuts5=cuts[5];

  string name3=folder+"/"+sort_file+"10.dat";
  ofstream PRINT3(name3.c_str());
	if (!PRINT3.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

  for (int c=0;c<99;c=c+10) PRINT3 << c << " " << c+10 << " " <<  cuts[c]  << " "  << cuts[c+10] << endl;
  double cuts10=cuts[10];
  PRINT3.close();

  //cout << cuts[1] << " " << cuts10 << " " << cuts1 << endl;

  return cuts10;

}

int sortmake2(string folder,std::vector<eccs> all){

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax,cmin;


  std::sort (all.begin(), all.end(), compareByNpart);

  int dbin=all.size()/100;

	vector<int> cuts;
	for (int c=0;c<100;c++) cuts.push_back(all[c*dbin].npart);

	int end=all.size()-1;
	cuts[100]=all[end].npart;



  string name=folder+"/npartcen1.dat";
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}

  for (int c=0;c<99;c++) PRINT << c << " " << c+1 << " " <<  cuts[c]  << " "  << cuts[c+1] << endl;
    double cuts1=cuts[1];
  PRINT.close();


  string name2=folder+"/npartcen5.dat";
  ofstream PRINT2(name2.c_str());
	if (!PRINT2.is_open())
 	{
 	cout << "Can't open " << name2 << endl;
 	exit(1);
 	}

  for (int c=0;c<99;c=c+5) PRINT2 << c << " " << c+5 << " " <<  cuts[c]  << " "  << cuts[c+5] << endl;
   double cuts5=cuts[5];
  PRINT2.close();

  string name3=folder+"/npartcen10.dat";
  ofstream PRINT3(name3.c_str());
	if (!PRINT3.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

  for (int c=0;c<99;c=c+10) PRINT3 << c << " " << c+10 << " " <<  cuts[c]  << " "  << cuts[c+10] << endl;
   double cuts10=cuts[10];
  PRINT3.close();

	//cout << cuts[1] << " " << cuts10 << " " << cuts1 << endl;

  return cuts10;

}

void sort(string folder, string sort_file, std::vector<eccs> all,  std::vector< std::vector<eccs> > & out, std::vector<double> & cens,int bin){

//      calculates max and min of each centrality window according to binning
	vector <double> clist;
	vector<double> cmax,cmin;


	for (int i=0;i<(100-bin);i+=bin) {

	cmax.push_back(i+bin);
	cmin.push_back(i);
	cens.push_back(i+bin/2.);
	}
	cout << "check sizes cmin " << cmin.size() << " " << cmax.size() << endl;

  int censize=cens.size();
	vector<double> nph,npl;
	nph.resize(censize);
	npl.resize(censize);


	string named=folder+"/"+sort_file+".dat";
	//cout << "named  " << named << endl;
	ifstream inputd(named.c_str());
	if (!inputd.is_open())
 	{
 	cout << "Can't open " << named << endl;
 	exit(1);
 	}
 	string line;
 	int imin=0,imax=0;

	while (getline(inputd,line)) {

	std::vector<double> y (2,0) ;
	std::vector<double> y2 (4,0) ;
	std::vector<std::string> x = split(line, ' ');
	if (x.size()<4) std::vector<std::string> x = split(line, '\t');

	//cout << "x.size  " << x.size() << endl;

	for(int j=0;j<2;j++)
	{
	stringstream s;
	s << x[j];
	s >> y[j];
	}
	for(int j=2;j<4;j++)
	{
	stringstream s;
	s << x[j];
	s >> y2[j];
	}

	//cout <<"imin  " << imin << " imax  " << imax  << " nph  " << nph.size() << " npl  " << npl.size() << endl;
	if (cmin[imin]==y[0]){//min of the centrality bin
	 nph[imin]=y2[2]; //max entropy for centrality bin
	 imin++;
	 if (imin>censize) break;
	 }
	if (cmax[imax]==y[1]){ //max of the centrality bin
	 npl[imax]=y2[3];//min entropy for centrality bin
	 imax++;
	 if (imax>censize) break;
	 }



	}
	inputd.close();



	//sets up variables to the approparite number of bins in centrality
	int tot=all.size();

	out.resize(censize);
	for (int i=0;i<tot;i++){
		for (int c=0;c<censize;c++){
			if (all[i].si<=nph[c]&&all[i].si>npl[c]){
			out[c].push_back(all[i]);
			break;
			}
		}
	}


}


void NPsc(string name2,  string chg,pars p,  std::vector<eccs>  & ec){
  // sort eccentricities by Npart, store in list

  int max=maxnp(ec);
  int size=ec.size();
  
  std::vector< std::vector<eccs> > list;
  list.resize(max+1);
  for (int i=0;i<size;i++){
  	list[ec[i].npart].push_back(ec[i]);
  }
  
  std::vector< double > sc32,sc42,sc43,v2,v3,v4;
  sc32.resize(max+1,0);
  sc42.resize(max+1,0);
  sc43.resize(max+1,0);
  v2.resize(max+1,0);
  v3.resize(max+1,0);
  v4.resize(max+1,0);
  
   for (int n=0;n<=max;n++){
  	int full=list[n].size();
  	for (int f=0;f<full;f++){
  		v2[n]+=pow(list[n][f].ec[2],2);
  		v3[n]+=pow(list[n][f].ec[3],2);
  		v4[n]+=pow(list[n][f].ec[4],2);
  		sc32[n]+=pow(list[n][f].ec[2]*list[n][f].ec[3],2);
  		sc42[n]+=pow(list[n][f].ec[2]*list[n][f].ec[4],2);
 		sc43[n]+=pow(list[n][f].ec[4]*list[n][f].ec[3],2);
  	}
  	v2[n]/=full;
  	v3[n]/=full;
  	v4[n]/=full;
  	sc32[n]/=full;
  	sc42[n]/=full;
  	sc43[n]/=full;
  }
  //cout << "list size 0 " << list[0].size() << endl << "list size 1 " << list[1].size() << endl << "list size max " << list[max].size() << endl << "list size " << list.size() << endl;
  string name=name2+"/obsb3_"+chg+"/eNSC_Npart.dat";
  cout << "max " << max << endl;
  ofstream fileout(name.c_str());
  cout << "Print is open" << endl;
	if (!fileout.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}
	//cout << "testing" << endl;
	for (int n=1;n<=max;n++){
	//cout << n << endl;
    fileout << n << " " << sc32[n]/(v2[n]*v3[n])-1  << " " << sc42[n]/(v2[n]*v4[n])-1  << " "<< sc43[n]/(v4[n]*v3[n])-1  <<  endl;
  }

  fileout.close();
	cout << "done np sc" << endl;
}


void ultracen(string name2,  string chg,std::vector<eccs> all, double scut, int npcut){
  int tot=all.size();
  std::vector<eccs> ssort,npsort;
  for (int i=0;i<tot;i++){
    if (all[i].si>=scut) ssort.push_back(all[i]);
    if (all[i].npart>=npcut) npsort.push_back(all[i]);
  }

cout << ssort.size() << " " << npsort.size() << endl;

  std::sort(ssort.begin(), ssort.end(), compareByLength);
  std::sort(npsort.begin(), npsort.end(), compareByLength);

 string sn="entropy";
  ucprint(name2,chg,sn,ssort);
   sn="npart";
  ucprint(name2,chg,sn,npsort);
  sn="entropy";
  ucprintCGC(name2,chg,sn,ssort);
  sn="npart";
  ucprintCGC(name2,chg,sn,npsort);
	cout << "mpart" << endl;
}

void ucprint(string name2,string chg,string pname,vector<eccs> ssort){
  int slen=ssort.size();

  int sbin=(slen/20);

  int sc=0;
  int npc=0;
  cout << "starting" << endl;
  vector<double> v2,v3,sc32,M;
  double mv2=0,mv3=0,msc32=0,mM=0;
  while(sc<slen){
    double c2=0,c3=0,c32=0,cM=0;
    for(int i=sc;i<(sc+sbin);i++){
      if (i>=slen) break;
      c2+=pow(ssort[i].v[2],2);
      c3+=pow(ssort[i].v[3],2);
      c32+=pow(ssort[i].v[2]*ssort[i].v[3],2);
      cM+=ssort[i].M;
    }
    mv2+=c2;
    mv3+=c3;
    msc32+=c32;
    mM+=cM;
    v2.push_back(sqrt(c2/sbin));
    v3.push_back(sqrt(c3/sbin));
    sc32.push_back(c32/sbin);
    M.push_back(cM/sbin);

  sc+=sbin;
  }

  double c2=0,c3=0,c32=0,cM=0;
  for(int i=sc;i<slen;i++){
    c2+=pow(ssort[i].v[2],2);
    c3+=pow(ssort[i].v[3],2);
    c32+=pow(ssort[i].v[2]*ssort[i].v[3],2);
    cM+=ssort[i].M;
  }
  mv2+=c2;
  mv3+=c3;
  msc32+=c32;
  mM+=cM;
  v2.push_back(sqrt(c2/sbin));
  v3.push_back(sqrt(c3/sbin));
  sc32.push_back(c32/sbin);
  M.push_back(cM/sbin);


  mv2/=slen;
  mv3/=slen;
  msc32/=slen*mv2*mv3;
  mM/=slen;

  cout << slen << endl;

  int olen=v2.size();
  vector<double> erv2,erv3,ersc32;
  for(int i=0;i<olen;i++){
    double sv2=0,sv3=0,ssc32=0;


    for(int j=i*sbin;j<((i+1)*sbin);j++){
      if (j>=slen) break;
      sv2+=pow(v2[i]/sqrt(mv2)-sqrt((pow(v2[i],2)*sbin-pow(ssort[j].v[2],2))/(sbin-1))/sqrt(mv2),2);
      sv3+=pow(v3[i]/sqrt(mv3)-sqrt((pow(v3[i],2)*sbin-pow(ssort[j].v[3],2))/(sbin-1))/sqrt(mv3),2);
      ssc32+=pow(sc32[i]/(pow(v2[i],2)*pow(v3[i],2))-(sbin-1)*(sc32[i]*sbin-pow(ssort[i].v[2]*ssort[i].v[3],2))/((pow(v2[i],2)*sbin-pow(ssort[j].v[2],2))*(pow(v3[i],2)*sbin-pow(ssort[j].v[3],2))),2);
    }
    erv2.push_back(sqrt((sbin-1)*sv2/sbin));
    erv3.push_back(sqrt((sbin-1)*sv2/sbin));
    ersc32.push_back(sqrt((sbin-1)*ssc32/sbin));
  }
  cout << "starting3" << endl;

  string name=name2+"/obsb3_"+chg+"/ultracen_"+pname+".dat";
  cout << name << endl;
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}


  for(int i=0;i<olen;i++){
    PRINT << M[i]/mM << " " << v2[i]/sqrt(mv2)  << " " << erv2[i] << " "<< v3[i]/sqrt(mv3) << " " << erv3[i] <<   " " <<  sc32[i]/(pow(v2[i],2)*pow(v3[i],2))-1 << " " << ersc32[i] << endl;
  }

  cout << "0-1% average multiplicity=" << mM << endl;

  PRINT.close();


}

void ucprintCGC(string name2,string chg,string pname,vector<eccs> ssort){
  int slen=ssort.size();

  int sbin=(slen/20);

  int sc=0;
  int npc=0;
  vector<double> v2,v3,sc32,M,v24v22,v34v32,cn2,cn3;
  double mv2=0,mv3=0,msc32=0,mM=0;

  double amean2=0,amean3=0;
  vector <double> norm,Inot;
  double I0=0;
  while(sc<slen){
    double c2=0,c3=0,c32=0,cM=0;
    double mean2=0,mean3=0;

    norm.push_back(0);
    Inot.push_back(0);
    int cur=norm.size()-1;
    double I0i=0;
    for(int i=sc;i<(sc+sbin);i++){
      if (i>=slen) break;
      c2+=pow(ssort[i].vND[2],2);
      c3+=pow(ssort[i].vND[3],2);
      c32+=ssort[i].vND[2]*ssort[i].vND[3];
      mean2+=ssort[i].vND[2];


      mean3+=ssort[i].vND[3];
      cM+=ssort[i].M;
      I0i+=ssort[i].M*ssort[i].M;
      norm[cur]++;

    }
    mv2+=mean2;
    mv3+=mean3;
    msc32+=c32;
    mM+=cM;
    amean2+=mean2;
    amean3+=mean3;
    v2.push_back(sqrt(mean2/pow(I0i,1)));
    v3.push_back(sqrt(mean3/pow(I0i,1)));
    sc32.push_back(c32/pow(I0i,2));
    M.push_back(cM/norm[cur]);
    cn2.push_back(c2);
    cn3.push_back(c3);
    Inot[cur]=I0i;
    I0+=I0i;

    double s2=-2.*(c2/pow(I0i,2)-pow(mean2/pow(I0i,1),1))/pow(mean2/pow(I0i,1),2);
    v24v22.push_back(s2);
    double s3=-2.*(c3/pow(I0i,2)-pow(mean3/pow(I0i,1),1))/pow(mean3/pow(I0i,1),2);
    v34v32.push_back(s3);
    cout << v2[cur] << " " <<norm[cur] << " "<< sqrt(mean2/pow(I0i,1)) << " " << sqrt(mean3/pow(I0i,2)) << " " << s2 << endl;
  sc+=norm[cur];
  }




  mv2/=pow(I0,1);
  mv3/=pow(I0,1);
  amean2/=pow(I0,1);
  amean3/=pow(I0,1);
  msc32/=pow(I0,2)*amean2*amean3;
  mM/=slen;

	cout << "averaged" << endl;


  int olen=v2.size();
  vector<double> erv2,erv3,ersc32,errv24v22,errv34v32;
  for(int i=0;i<olen;i++){
    double sv2=0,sv3=0,ssc32=0;
    double vs[4]={0},v4s[4]={0},ms[4]={0};
    for(int j=i*sbin;j<((i+1)*sbin);j++){
      if (j>=slen) break;
      double Msub=Inot[i]-pow(ssort[j].M,1);
      sv2+=pow(v2[i]/sqrt(mv2)-sqrt((pow(v2[i],2)*pow(Inot[i],1)-ssort[j].vND[2])/pow(Msub,1))/sqrt(mv2),2);
      sv3+=pow(v3[i]/sqrt(mv3)-sqrt((pow(v3[i],2)*pow(Inot[i],1)-ssort[j].vND[3])/pow(Msub,1))/sqrt(mv3),2);
      ssc32+=pow(sc32[i]/(pow(v2[i],2)*pow(v3[i],2))-(sc32[i]*pow(Inot[i],2)-ssort[i].vND[2]*ssort[i].vND[3])/((pow(v2[i],2)*pow(Inot[i],1)-ssort[j].vND[2])*(pow(v3[i],2)*pow(Inot[i],1)-ssort[j].v[3])),2);



       vs[2]=cn2[i]-ssort[j].vND[2]*ssort[j].vND[2];
       ms[2]=pow(v2[i],2)*pow(Inot[i],1)-ssort[j].vND[2];
       vs[3]=cn3[i]-ssort[j].vND[3]*ssort[j].vND[3];
       ms[3]=pow(v3[i],2)*pow(Inot[i],1)-ssort[j].vND[3];

       double sub=-2.*(vs[2]/pow(Msub,1)-pow(ms[2]/pow(Msub,1),2))/(pow(ms[2]/pow(Msub,1),2));
       v4s[2]+=pow(v24v22[i]-sub,2);
       sub=-2.*(vs[3]/pow(Msub,1)-pow(ms[3]/pow(Msub,2),2))/(pow(ms[3]/pow(Msub,1),2));
       v4s[3]+=pow(v34v32[i]-sub,2);
    }
    erv2.push_back(sqrt((norm[i]-1)*sv2/norm[i]));
    erv3.push_back(sqrt((norm[i]-1)*sv2/norm[i]));
    ersc32.push_back(sqrt((norm[i]-1)*ssc32/norm[i]));
    errv24v22.push_back(sqrt((norm[i]-1)*v4s[2]/norm[i]));
    errv34v32.push_back(sqrt((norm[i]-1)*v4s[3]/norm[i]));
  }

  cout << "done" << endl;

  string name=name2+"/obsb3_"+chg+"/ultracen_"+pname+"CGC.dat";
  ofstream PRINT(name.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}


  for(int i=0;i<olen;i++){
    PRINT << M[i]/mM << " " << v2[i]/sqrt(mv2)  << " " << erv2[i] << " "<< v3[i]/sqrt(mv3) << " " << erv3[i] <<   " " <<  sc32[i]/(pow(v2[i],2)*pow(v3[i],2))-1 << " " << ersc32[i] << " "  << v24v22[i] << " " << errv24v22[i] << " "   << v34v32[i] << " "  << errv34v32[i] <<  endl;
  }

  PRINT.close();
  cout << "printed" << endl;

}


//
//
// void sort(string folder, string sort_file, std::vector<eccs> all,  std::vector< std::vector<eccs> > & out, std::vector<double> & cens,int bin){
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
// 	std::vector<double> y (2,0) ;
// 	std::vector<double> y2 (4,0) ;
// 	std::vector<std::string> x = split(line, ' ');
// 	if (x.size()<4) std::vector<std::string> x = split(line, '\t');
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


double pred(double k1, double k2, double ec){
	return (k1*ec+k2*pow(ec,3));
}


void vns(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin){

	//string name3="trento/"+name2+"/ND/vns1.dat";
	string name3=name2+"/obsb3_"+chg+"/vns1bin"+convertInt(bin)+"rat.dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

  string namet4=name2+"/obsb3_"+chg+"/vns1bin"+convertInt(bin)+"to4.dat";
	ofstream PRINTt4(namet4.c_str());
	if (!PRINTt4.is_open())
 	{
 	cout << "Can't open " << namet4 << endl;
 	exit(1);
 	}

	int ssize=sall.size();

	for (int j=0;j<ssize;j++ ){
		int evs=sall[j].size();
		double v[7]={0},v4[7]={0},v4t4[7]={0},v6[7]={0},sv4[7]={0},sv6[7]={0},err[7]={0},errt4[7]={0},err6[7]={0};
		for(int ev=0;ev<evs;ev++){
		 for (int n=2;n<=5;n++) {
		 	v[n]+=sall[j][ev].ec[n]*sall[j][ev].ec[n];
		 	sv4[n]+=pow(sall[j][ev].ec[n],4);
		 	sv6[n]+=pow(sall[j][ev].ec[n],6);
		 }
		}
		for (int n=2;n<=5;n++) {
		v4[n]=pow(2.*v[n]/evs*v[n]/evs-sv4[n]/evs,0.25)/sqrt(v[n]/evs);
    v4t4[n]=(2.*v[n]/evs*v[n]/evs-sv4[n]/evs)/pow(v[n]/evs,2);
		v6[n]=pow(0.25*(sv6[n]/evs-9.*v[n]/evs*sv4[n]/evs+12.*pow(v[n]/evs,3) ),1./6.)/pow(2.*v[n]/evs*v[n]/evs-sv4[n]/evs,0.25);
		}


		double vs[7]={0},v4s[7]={0},v4st4[7]={0},sv4s[7]={0},v6s[7]={0},sv6s[7]={0};
		for(int ev=0;ev<evs;ev++){
		 for (int n=2;n<=5;n++) {
		 	vs[n]=v[n]-sall[j][ev].ec[n]*sall[j][ev].ec[n];
		 	sv4s[n]=sv4[n]-pow(sall[j][ev].ec[n],4);
		 	sv6s[n]=sv6[n]-pow(sall[j][ev].ec[n],6);
		 }

		 for (int n=2;n<=5;n++) {
		 	v4s[n]+=pow(v4[n]-pow(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1),0.25)/sqrt(vs[n]/(evs-1)),2);
      v4st4[n]+=pow(v4t4[n]-(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1))/pow(vs[n]/(evs-1),2),2);
			v6s[n]+=pow(v6[n]-pow(0.25*(sv6s[n]/(evs-1)-9.*vs[n]/(evs-1)*sv4s[n]/(evs-1)+12.*pow(vs[n]/(evs-1),3) ),1./6.)/pow(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1),0.25),2);


			}

		}
		for (int n=2;n<=5;n++) {
			err[n]=sqrt((evs-1)*v4s[n]/evs);
      errt4[n]=sqrt((evs-1)*v4st4[n]/evs);
			err6[n]=sqrt((evs-1)*v6s[n]/evs);
			}

		PRINT << cens[j] << " " ;

		for (int n=2;n<=5;n++) PRINT <<  v4[n] << " " << err[n] << " ";
		for (int n=2;n<=5;n++) PRINT <<  v6[n] << " " << err6[n] << " ";
		PRINT << endl;

    PRINTt4 << cens[j] << " " ;

    for (int n=2;n<=5;n++) PRINTt4 <<  v4t4[n] << " " << errt4[n] << " ";
    PRINTt4 << endl;


	}

	PRINT.close();
  PRINTt4.close();
}


void radii(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin){

	
	string name3=name2+"/obsb3_"+chg+"/radii_"+convertInt(bin)+".dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

	int ssize=sall.size();

	for (int j=0;j<ssize;j++ ){
		int evs=sall[j].size();
		double r2=0;
		for(int ev=0;ev<evs;ev++){
		  r2+=sall[j][ev].R*sall[j][ev].R;
		}
		r2/=evs;
		PRINT << cens[j] << " " << sqrt(r2) << endl;
	}

	PRINT.close();
}

void vCGC(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin){

	//string name3="trento/"+name2+"/ND/vns1.dat";
	string name3=name2+"/obsb3_"+chg+"/CGCfluc"+convertInt(bin)+".dat";
	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

	int ssize=sall.size();

	double v2cen=0,v3cen=0,Icen=0;
	for (int j=0;j<ssize;j++ ){
		int evs=sall[j].size();
		double v[4]={0},v4[4]={0},vnm=0,err[4]={0},mv[4]={0},v2[4]={0};
		double I0i=0;
		for(int ev=0;ev<evs;ev++){
		 I0i+=sall[j][ev].M*sall[j][ev].M;
		 for (int n=2;n<=3;n++) {
      			mv[n]+=sall[j][ev].vND[n];
		 	v[n]+=sall[j][ev].vND[n]*sall[j][ev].vND[n];

		 }
    		 vnm+=sall[j][ev].vND[2]*sall[j][ev].vND[3];
		}
		for (int n=2;n<=3;n++) v4[n]=-2.*(v[n]/evs-pow(mv[n]/evs,2))/pow(mv[n]/evs,2);

    vnm/=evs;
    double sc32=vnm/(mv[2]/evs*mv[3]/evs)-1;
    for (int n=2;n<=3;n++) v2[n]=sqrt(mv[n]/I0i);
    if (j==0) {
    	v2cen=v2[2];
    	v3cen=v2[3];
    	Icen= I0i;
    }


		double vs[4]={0},v4s[4]={0},ms[4]={0};
    double ssc32=0;
		for(int ev=0;ev<evs;ev++){
		 for (int n=2;n<=3;n++) {
		 	vs[n]=v[n]-sall[j][ev].vND[n]*sall[j][ev].vND[n];
               ms[n]=mv[n]-sall[j][ev].vND[n];
		 }

		 for (int n=2;n<=3;n++) {
       double sub=-2.*(vs[n]/(evs-1)-pow(ms[n]/(evs-1),2))/(pow(ms[n]/(evs-1),2));
		 	v4s[n]+=pow(v4[n]-sub,2);
			}
      ssc32+=pow(sc32+1-(evs-1)*(vnm*evs-sall[j][ev].vND[2]*sall[j][ev].vND[3])/((vs[2]*evs-sall[j][ev].vND[2])*(vs[3]*evs-sall[j][ev].vND[3])),2);

		}
		for (int n=2;n<=3;n++) {
			err[n]=sqrt((evs-1)*v4s[n]/evs);

			}
      err[0]=sqrt((evs-1)*ssc32/evs);
		PRINT << cens[j] << " " << I0i/Icen << " " <<  v2[2]/v2cen <<  " "  <<  v2[3]/v3cen <<  " ";
		for (int n=2;n<=3;n++) PRINT <<  v4[n] << " " << err[n] << " ";
		PRINT << sc32 << " " << err[0] <<  endl;


	}

	PRINT.close();
}
//void vnsb1(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens){

//	//string name3="trento/"+name2+"/ND/vns1.dat";
//	string name3="trento/"+name2+"/obsb3_"+chg+"/vns1bin1n.dat";
//	ofstream PRINT(name3.c_str());
//	if (!PRINT.is_open())
// 	{
// 	cout << "Can't open " << name3 << endl;
// 	exit(1);
// 	}

//	int ssize=sall.size();
//
//	for (int j=0;j<ssize;j++ ){
//		int evs=sall[j].size();
//		double v[7]={0},v4[7]={0},v6[7]={0},sv4[7]={0},sv6[7]={0},err[7]={0},err6[7]={0};
//		for(int ev=0;ev<evs;ev++){
//		 for (int n=2;n<=5;n++) {
//		 	v[n]+=sall[j][ev].ec[n]*sall[j][ev].ec[n];
//		 	sv4[n]+=pow(sall[j][ev].ec[n],4);
//		 	sv6[n]+=pow(sall[j][ev].ec[n],6);
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
//		 	vs[n]=v[n]-sall[j][ev].ec[n]*sall[j][ev].ec[n];
//		 	sv4s[n]=sv4[n]-pow(sall[j][ev].ec[n],4);
//		 	sv6s[n]=sv6[n]-pow(sall[j][ev].ec[n],6);
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



void vnsnorat(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens,int bin){

//	string name3="trento/"+name2+"/ND/vns1bin5.dat";
//	string name4="trento/"+name2+"/ND/vn2err.dat";
	string name3=name2+"/obsb3_"+chg+"/enstotbin"+convertInt(bin)+".dat";
	string name4=name2+"/obsb3_"+chg+"/en2errbin"+convertInt(bin)+".dat";
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

	int ssize=sall.size();

	for (int j=0;j<ssize;j++ ){
		int evs=sall[j].size();
		double v[7]={0},v2[7]={0},v4[7]={0},v6[7]={0},sv4[7]={0},sv6[7]={0},err2[7]={0},err[7]={0},err6[7]={0};
		double mpt=0, mpt2=0,rpt=0;
		double v2v3=0;

		for(int ev=0;ev<evs;ev++) mpt+=sall[j][ev].pT;
		mpt/=evs;

		for(int ev=0;ev<evs;ev++){
		 for (int n=2;n<=5;n++) {
		 	v[n]+=pow(sall[j][ev].ec[n],2);
		 	sv4[n]+=pow(sall[j][ev].ec[n],4);
		 	sv6[n]+=pow(sall[j][ev].ec[n],6);
		 }

		 mpt2+=pow(sall[j][ev].pT-mpt,2.);
		}
		for (int n=2;n<=5;n++) {
		v2[n]=sqrt(v[n]/evs);
		v4[n]=2.*v[n]/evs*v[n]/evs-sv4[n]/evs;
		v6[n]=0.25*(sv6[n]/evs-9.*v[n]/evs*sv4[n]/evs+12.*pow(v[n]/evs,3) );
		}
		v2v3=sqrt(v[2]/v[3]);

		mpt2/=evs;
		rpt=sqrt(mpt2)/mpt;

		double vs[7]={0},v2s[7]={0},v4s[7]={0},sv4s[7]={0},v6s[7]={0},sv6s[7]={0},sv2v3=0,ev2v3=0;
		double empt=0, empt2=0,erpt=0;
		double smpt=0, smpt2=0,srpt=0;
		for(int ev=0;ev<evs;ev++){
		 for (int n=1;n<=6;n++) {
		 	vs[n]=v[n]-pow(sall[j][ev].ec[n],2);
		 	sv4s[n]=sv4[n]-pow(sall[j][ev].ec[n],4);
		 	sv6s[n]=sv6[n]-pow(sall[j][ev].ec[n],6);
		 }

		 for (int n=2;n<=5;n++) {
		 	v2s[n]+=pow(v2[n]-sqrt(vs[n]/(evs-1)),2);
			v4s[n]+=pow(v4[n]-(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1) ),2);
			v6s[n]+=pow(v6[n]-0.25*(sv6s[n]/(evs-1)-9.*vs[n]/(evs-1)*sv4s[n]/(evs-1)+12.*pow(vs[n]/(evs-1),3) ),2);
			}
			double subpt=(mpt*evs-sall[j][ev].pT)/(evs-1.);

			double subpt2=(mpt2*evs-pow(sall[j][ev].pT-subpt,2) )/(evs-1);

			sv2v3+=pow(v2v3-sqrt(vs[2]/vs[3]),2);
		 	smpt+=pow(mpt-subpt,2);
		 	smpt2+=pow(mpt2-subpt2,2);
		 	srpt=pow(rpt-sqrt(subpt2)/subpt,2);
		}



		for (int n=2;n<=5;n++) {
			err2[n]=sqrt((evs-1)*v2s[n]/evs);
			err[n]=sqrt((evs-1)*v4s[n]/evs);
			err6[n]=sqrt((evs-1)*v6s[n]/evs);
		}
		ev2v3=sqrt((evs-1)*sv2v3/evs);
		empt=sqrt((evs-1)*smpt/evs);
		empt2=sqrt((evs-1)*smpt2/evs);
		erpt=sqrt((evs-1)*srpt/evs);

		PRINT << cens[j] << " " ;

		for (int n=2;n<=5;n++) PRINT <<  v4[n] << " " << err[n] << " ";
		for (int n=2;n<=5;n++) PRINT <<  v6[n] << " " << err6[n] << " ";
		PRINT << v2v3 << " " << ev2v3 << endl;

		PRINT2 << cens[j] << " " ;
		for (int n=2;n<=5;n++) PRINT2 <<  v2[n] << " " << err2[n] << " ";
		PRINT2 << endl;
	}

	PRINT.close();
	PRINT2.close();
}


void predcum(string full, string runtype,  string runname,string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens, int bin){


// Reads in the previously calculated kappa's for the linear+cupic response
    string name=runtype+"/obschgv1/nonlin_bin"+convertInt(bin)+".dat";
    ifstream input(name.c_str());
    if (!input.is_open())
      {
      cout << "Can't open predictions file " << name << endl;
      return;
      }

    string line;
    int cc=0;
    vector <double> k1v2,k2v2;
    vector <double> k1v3,k2v3;
    while (getline(input,line)) {


    std::vector<double> y (13,0) ;
    std::vector<std::string> x = split(line, ' ');


    for(int j=0;j<13;j++)
    {
    stringstream s;
    s << x[j];
    s >> y[j];
    }
    k1v2.push_back(y[4-1]);
    k2v2.push_back(y[6-1]);
    k1v3.push_back(y[10-1]);
    k2v3.push_back(y[12-1]);
    cc++;
    }
    input.close();
// Done reading in kappas



// define names of outputfiles
  string name6=full+"/obsb3_"+chg+"/predictv24v22_bin"+convertInt(bin)+"_"+runname+".dat";

  ofstream PRINT4(name6.c_str());
  if (!PRINT4.is_open())
  {
  cout << "Can't open " << name6 << endl;
  exit(1);
  }

  string name2=full+"/obsb3_"+chg+"/predictv22_bin"+convertInt(bin)+"_"+runname+".dat";

  ofstream PRINT(name2.c_str());
  if (!PRINT.is_open())
  {
  cout << "Can't open " << name2 << endl;
  exit(1);
  }
  //end creating outputfiles

//calculates vn{2} and vn{4}/vn{2} using linear+cubic response from the eccentricities (no actual flow used here)
// j is the centrality bin
int ssize=sall.size();
	for (int j=0;j<ssize;j++ ){
		int evs=sall[j].size();
		double err[7]={0}, err2[7]={0};

		//predicts v2{4}/v2{2} from eccentricities
		double ecc22[4]={0},ecc44[4]={0},pecv2[4]={0},pecv24[4]={0};
		for(int ev=0;ev<evs;ev++){
		 for (int n=2;n<=3;n++) {
      double en=0;
		 	if (n==2)en=pred(k1v2[j],k2v2[j],sall[j][ev].ec[n]);
      if (n==3)en=pred(k1v3[j],k2v3[j],sall[j][ev].ec[n]);
		 	ecc22[n]+=pow(en,2);
		 	ecc44[n]+=pow(en,4);
		 }// calculates moments of the distributions
		}
		for (int n=2;n<=3;n++) {
			ecc22[n]/=evs;
			ecc44[n]/=evs;
      pecv2[n]=sqrt(ecc22[n]);
			pecv24[n]=pow(2*pow(ecc22[n],2)-ecc44[n],0.25)/sqrt(ecc22[n]);
		} // these are the central values, error bars calculated below


   //here we used jackknife error to estimate the error bands
		double vs[7]={0},v2s[7]={0},v4s[7]={0},sv4s[7]={0};
		for(int ev=0;ev<evs;ev++){
		 for (int n=2;n<=3;n++) {
      double en=0;
      if (n==2)en=pred(k1v2[j],k2v2[j],sall[j][ev].ec[n]);
      if (n==3)en=pred(k1v3[j],k2v3[j],sall[j][ev].ec[n]);
		 	vs[n]=(ecc22[n]*evs-pow(en,2))/(evs-1);
		 	sv4s[n]=(ecc22[n]*evs-pow(en,4))/(evs-1);
		 }

		 for (int n=2;n<=3;n++) {
		 	v2s[n]+=pow(pecv2[n]-sqrt(vs[n]),2);
			v4s[n]+=pow(pecv24[n]-pow(2.*vs[n]*vs[n]-sv4s[n],0.25) ,2);
    }


		}


		for (int n=2;n<=3;n++) {
			err[n]=sqrt((evs-1)*v2s[n]/evs);
			err2[n]=sqrt((evs-1)*v4s[n]/evs);
		}
    //done calculating errors


    //printing off vn{4} and vm{2} in separate files


		PRINT4 << cens[j] << " " ;
		for (int n=2;n<=3;n++) PRINT4 << pecv24[n]  << " " << err2[n] << " " ;
		PRINT4 << endl;


    PRINT << cens[j] << " " ;
    for (int n=2;n<=3;n++) PRINT << pecv2[n]  << " " << err[n] << " " ;
    PRINT << endl;
    }


    PRINT.close();
    PRINT4.close();


}


//void vns3(string name2, string chg, std::vector< std::vector<eccs> > & sall, std::vector<double> & cens){

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

//	int ssize=sall.size();
//
//	for (int j=0;j<ssize;j++ ){
//		int evs=sall[j].size();
//		double v[7]={0},v2[7]={0},v4[7]={0},v6[7]={0},sv4[7]={0},sv6[7]={0},err2[7]={0},err[7]={0},err6[7]={0};
//		double mpt=0, mpt2=0,rpt=0;
//		double v2v3=0;
//
//		for(int ev=0;ev<evs;ev++) mpt+=sall[j][ev].pT;
//		mpt/=evs;
//
//		for(int ev=0;ev<evs;ev++){
//		 for (int n=2;n<=5;n++) {
//		 	v[n]+=pow(sall[j][ev].ec[n],2);
//		 	sv4[n]+=pow(sall[j][ev].ec[n],4);
//		 	sv6[n]+=pow(sall[j][ev].ec[n],6);
//		 }
//
//		 mpt2+=pow(sall[j][ev].pT-mpt,2.);
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
//		 	vs[n]=v[n]-pow(sall[j][ev].ec[n],2);
//		 	sv4s[n]=sv4[n]-pow(sall[j][ev].ec[n],4);
//		 	sv6s[n]=sv6[n]-pow(sall[j][ev].ec[n],6);
//		 }
//
//		 for (int n=2;n<=5;n++) {
//		 	v2s[n]+=pow(v2[n]-sqrt(vs[n]/(evs-1)),2);
//			v4s[n]+=pow(v4[n]-(2.*vs[n]/(evs-1)*vs[n]/(evs-1)-sv4s[n]/(evs-1) ),2);
//			v6s[n]+=pow(v6[n]-0.25*(sv6s[n]/(evs-1)-9.*vs[n]/(evs-1)*sv4s[n]/(evs-1)+12.*pow(vs[n]/(evs-1),3) ),2);
//			}
//			double subpt=(mpt*evs-sall[j][ev].pT)/(evs-1.);
//
//			double subpt2=(mpt2*evs-pow(sall[j][ev].pT-subpt,2) )/(evs-1);
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

