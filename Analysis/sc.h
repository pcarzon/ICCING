#ifndef SC_H_
#define SC_H_

#include <string>
#include <iostream>
#include <string.h>
#include <vector>

# include "global.h"

using namespace std;

class sc {
private:

	double fixsc[10],fixM4[10];
	double fixv6[10],fixv4[10],fixM6[10];


public:

	void eccsout(string name, vector< vector<eccs> > sall);

	double  recombM(std::vector< std::vector<eccs> > fix2, pars P1);

	double  recombMfix(vector< vector<eccs> > fix2, pars P1,int j) ;

	serr jack (  vector< vector<eccs> > sall, pars P1 );

	SCsub SCM(vector<eccs> vec, pars p1);

	SCsub EP3M3(vector<eccs> vec, pars p1);

	SCsub v4cum(vector<eccs> vec, pars p1);


	double EP3M(vector<eccs> vec, pars p1);

	double  recombMEP(std::vector< std::vector<eccs> > fix2, pars P1);

	double  recombMfixEP(vector< vector<eccs> > fix2, pars P1,int j) ;

	// jack knife event plane
	serr jackEP (  vector< vector<eccs> > sall, pars P1 );

	void runEP(string name, pars p1, vector< vector<eccs> > sall,int j,  std::vector<cent> cens);

	void run(string name, pars p1, vector< vector<eccs> > sall);
	void run(string name, pars p1, vector< vector<eccs> > sall,std::vector<cent> cens);

	void run(string name, pars p1, vector< vector<eccs> > sall,int cent);

	void runec(string name, pars p1, vector< vector<eccs> > sall);

	void scnobin (string name,   pars p1, vector< vector<eccs> > sall);

	double v2M(vector<eccs> vec, pars p1);

	double v2(vector<eccs> vec, pars p1);

	double v23M(vector<eccs> vec, pars p1);

	double v23(vector<eccs> vec, pars p1);

	void print(string oname,  vector<serr> out);

	void print(string oname,  vector<serr> out,std::vector<cent> cens);

	void print(string oname,  vector<serr> out,int j, int & first,  std::vector<cent> cens);
	void print(string oname,  vector<serr> out,int j, int & first);

	void print2(string oname,  vector<serr> out,int j, int & first);

	void correct(vector<eccs> & sall);
	void correct2(vector<eccs> & sall);
	void correctpb(vector<eccs> & sall, int j);

	serr  recombMfixvn(vector< vector<eccs> > fix2, pars P1,int j) ;
	serr  recombMvn(vector< vector<eccs> > fix2, pars P1);
	serr jackvn (  vector< vector<eccs> > sall, pars P1 );
	void runvn(string name, pars p1, vector< vector<eccs> > sall,int j);



	sc() {};
	~sc() {};
};



#endif
