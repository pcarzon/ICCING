#ifndef _SC_CPP_
#define _SC_CPP_

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
#include <iomanip>      // std::setprecision

# include "sc.h"
# include "vars.h"


using namespace std;

void sc::correctpb(vector<eccs> & sall, int j){

	int len=sall.size();
	int i=0;
	double mean=0;
	int k=0;
	while(i<len){


		//if (c==8) cout << sall[c][i].v[2] << endl;

		if (k<50000){
		if (sall[i].v[2]<=0||sall[i].v[2]>0.3||sall[i].M>10000||sall[i].M<0) {
		//cout << "v2 " << i << " " <<  sall[i].v[2] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[3]<=0||sall[i].v[3]>0.3) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
//		else if (sall[i].v[4]<=0||sall[i].v[4]>0.3) {
//		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
//		sall.erase( sall.begin() +i );
//		len--;
//		}
//		else if (sall[i].v[6]<=0||sall[i].v[6]>0.07) {
//		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
//		sall.erase( sall.begin() +i );
//		len--;
//		}
		else {
			mean+=sall[i].M;
			i++;
		}
		}
		else{
			if (sall[i].M>0||sall[i].v[2]<=0){
				mean+=sall[i].M;
				i++;
			}
			else {
				sall.erase( sall.begin() +i );
				len--;
			}
		}


		k++;
	}

	mean/=i;


	len=sall.size();
	i=0;
	while(i<len){
		if (j<55 && (sall[i].M>(1.2*mean)||sall[i].M<(0.8*mean) ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (j<60 && (sall[i].M>(1.3*mean)||sall[i].M<(0.7*mean) ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (j<65 && (sall[i].M>(1.4*mean)||sall[i].M<(0.6*mean) ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (j<80 && (sall[i].M>(1.75*mean)||sall[i].M<(0.25*mean) ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (j<95 && (sall[i].M>(1.9*mean)||sall[i].M<(0.15*mean) ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (j<=100 && sall[i].M<(0.20*mean ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else {
			i++;
		}

	}


}





void sc::correct(vector<eccs> & sall){

	int len=sall.size();
	int i=0;
	double mean=0;
	while(i<len){


		//if (c==8) cout << sall[c][i].v[2] << endl;

		if (sall[i].v[2]<=0||sall[i].v[2]>0.75||sall[i].M>10000||sall[i].M<0) {
		//cout << "v2 " << i << " " <<  sall[i].v[2] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[3]<=0||sall[i].v[3]>0.75) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[4]<=0||sall[i].v[4]>0.75) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[6]<=0||sall[i].v[6]>0.5) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else {
			mean+=sall[i].M;
			i++;
		}
	}

	mean/=i;

	len=sall.size();
	i=0;
	while(i<len){
		if (mean>50 && (sall[i].M>(1000*mean)||sall[i].M<(0*mean) ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else {
			i++;
		}

	}


}

void sc::correct2(vector<eccs> & sall){

	int len=sall.size();
	int i=0;
	double mean=0;
	while(i<len){


		//if (c==8) cout << sall[c][i].v[2] << endl;

		if (sall[i].v[2]<=0||sall[i].v[2]>0.7||sall[i].M>10000||sall[i].M<0) {
		//cout << "v2 " << i << " " <<  sall[i].v[2] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[3]<=0||sall[i].v[3]>0.7) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[4]<=0||sall[i].v[4]>0.7) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else if (sall[i].v[6]<=0||sall[i].v[6]>0.07) {
		//cout << "v3 " << i << " " <<  sall[i].v[3] << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else {
			mean+=sall[i].M;
			i++;
		}
	}

	mean/=i;


	len=sall.size();
	i=0;
	while(i<len){
		if (mean>0.2 && (sall[i].M>(1.1*mean)||sall[i].M<(0.9*mean)  ||sall[i].v[2]>0.4 ) ) {
		//cout << mean <<" " <<  sall[i].M << endl;
		sall.erase( sall.begin() +i );
		len--;
		}
		else {
			i++;
		}

	}


}



void sc::print(string oname,  vector<serr> out){

	string name3=oname+".dat";

	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

	int osize=out.size();
	for (int j=0;j<osize;j++) PRINT << out[j].cen << " " << out[j].ans << " " << out[j].err << endl;

	PRINT.close();
}

void sc::print(string oname,  vector<serr> out,std::vector<cent> cens){

	string name3=oname+".dat";

	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

	int osize=out.size();
	 PRINT << "Cent%\t dNdy\t Npart\t Observable\t err\n";
	for (int j=0;j<osize;j++) PRINT << cens[j].cls << "\t" << cens[j].M << "\t" <<cens[j].npart << "\t" << out[j].ans << " " << out[j].err << endl;

	PRINT.close();
}

void sc::print(string oname,  vector<serr> out, int j, int & first,  std::vector<cent> cens){
	string name3=oname+".dat";




	if (first==1) {

		ofstream PRINT(name3.c_str());
		if (!PRINT.is_open())
	 	{
	 	cout << "Can't open " << name3 << endl;
	 	exit(1);
	 	}
		 PRINT << "Cent%\t dNdy\t Npart\t Observable\t err\n";
		PRINT <<  cens[0].cls << "\t" << cens[0].M << "\t" <<cens[0].npart << "\t"<< out[0].ans << " " << out[0].err << endl;

		PRINT.close();
	}
	else {
		ofstream PRINT(name3.c_str(), std::ofstream::out | std::ofstream::app);
		if (!PRINT.is_open())
	 	{
	 	cout << "Can't open " << name3 << endl;
	 	exit(1);
	 	}
		int cs=j/10-1;
		PRINT <<  cens[cs].cls << "\t" << cens[cs].M << "\t" <<cens[cs].npart << "\t"<< out[0].ans << " " << out[0].err << endl;

		PRINT.close();
	}
}

void sc::print(string oname,  vector<serr> out, int j, int & first){
	string name3=oname+".dat";



	if (first==1) {

		ofstream PRINT(name3.c_str());
		if (!PRINT.is_open())
	 	{
	 	cout << "Can't open " << name3 << endl;
	 	exit(1);
	 	}

		PRINT << j << "\t"<< out[0].ans << " " << out[0].err << endl;

		PRINT.close();
	}
	else {
		ofstream PRINT(name3.c_str(), std::ofstream::out | std::ofstream::app);
		if (!PRINT.is_open())
	 	{
	 	cout << "Can't open " << name3 << endl;
	 	exit(1);
	 	}

		PRINT <<  j << "\t"<< out[0].ans << " " << out[0].err << endl;

		PRINT.close();
	}
}

void sc::print2(string oname,  vector<serr> out, int j, int & first){
	string name3=oname+".dat";

	if (first==1) {

		ofstream PRINT(name3.c_str());
		if (!PRINT.is_open())
	 	{
	 	cout << "Can't open " << name3 << endl;
	 	exit(1);
	 	}

		PRINT << out[0].cen << " " << out[0].ans << " " << out[0].err << " " << out[0].ans2 << " " << out[0].err2 << endl;

		PRINT.close();
	}
	else {
		ofstream PRINT(name3.c_str(), std::ofstream::out | std::ofstream::app);
		if (!PRINT.is_open())
	 	{
	 	cout << "Can't open " << name3 << endl;
	 	exit(1);
	 	}

		PRINT << out[0].cen << " " << out[0].ans << " " << out[0].err << " " << out[0].ans2 << " " << out[0].err2 << endl;

		PRINT.close();
	}
}

void sc::eccsout(string name, vector< vector<eccs> > sall){


	string name3=name+"eccs.dat";

	ofstream PRINT(name3.c_str());
	if (!PRINT.is_open())
 	{
 	cout << "Can't open " << name3 << endl;
 	exit(1);
 	}

	int asz=sall.size();
	for (int j=0;j<asz;j++){

	double ec2=0,ec3=0;
	int sz=sall[j].size();
	for (int i=0;i<=sz;i++){
		ec2+=sall[j][i].v[2];
		ec3+=sall[j][i].v[3];
	}
		ec2/=sz;
		ec3/=sz;

		PRINT << j << " " << ec2 << " " << ec3 << endl;

	}


	PRINT.close();


}


void sc::scnobin (string name,   pars P1, vector< vector<eccs> > sall){

	int end=sall.size();

	vector<serr> out;
	out.resize(end);



	double m=0,tot=0;
	for(int j=0;j<end;j++){
	SCsub ssub=sc::SCM(sall[j],P1);
	out[j].ans=ssub.vnvm;
	cout << out[j].ans << endl;

	int len=sall[j].size();
	tot+=len;
	for(int i=0;i<len;i++){
	  vector<eccs>  sub=sall[j];
	  sub.erase( sub.begin() +i );
	  SCsub ssub2=sc::SCM(sub,P1);
	  double fsub=ssub2.vnvm;
	  m+=(fsub-out[j].ans)*(fsub-out[j].ans);

	}
	out[j].err=sqrt((tot-1)*m/tot);
	}



 	print(name,out);

}


void sc::run(string name, pars p1, vector< vector<eccs> > sall,  std::vector<cent> cens){


	int end=sall.size();
	vector<serr> sc32all;
	int csub=0;
	for (int j=10;j<=end;j=j+10){

	vector< vector<eccs> > sub(&sall[j-10],&sall[j]);

	sc32all.push_back(sc::jack (sub,p1) );
	sc32all[csub].cen=j-5;
	//cout << sc32all[csub].cen << " " << sc32all[csub].ans << " " << sc32all[csub].err << endl;
	csub++;
	}

	print(name,sc32all,cens);

}

void sc::run(string name, pars p1, vector< vector<eccs> > sall,int j){



	vector<serr> sc32all;
	int csub=0;
	//for (int j=10;j<=60;j=j+10){

	vector< vector<eccs> > sub(&sall[j-10],&sall[j]);
	sc32all.push_back(sc::jack (sub,p1) );
	sc32all[csub].cen=j-5;
//	cout << sc32all[csub].cen << " " << sc32all[csub].ans << " " << sc32all[csub].err << endl;
//	getchar();
	csub++;
	//}

	int first=1;
	if (j>10) first=0;
	print(name,sc32all,j,first);

}

void sc::runec(string name, pars p1, vector< vector<eccs> > sall){

	int end=sall.size();
	vector<serr> sc32all;
	sc32all.resize(end);
	int csub=0;

	double vnvm=0;
	double totm4=0;
	for (int j=0;j<end;j++){
	SCsub fix=SCM(sall[j], p1);
	sc32all[j].ans=fix.vnvm;
	sc32all[j].err=0;
	sc32all[j].cen=j;
	}

	print(name,sc32all);
}



serr sc::jack (  vector< vector<eccs> > sall, pars P1 ){


	serr out;
	out.ans=sc::recombM(sall,P1);
	double m=0,tot=0;
	for(int j=0;j<10;j++){
	int len=sall[j].size();
	tot+=len;
	for(int i=0;i<len;i++){
	  vector< vector<eccs> > sub=sall;
	  sub[j].erase( sub[j].begin() +i );
	  double fsub=sc::recombMfix(sub,P1,j);
	  m+=(fsub-out.ans)*(fsub-out.ans);

	}}
	out.err=sqrt((tot-1)*m/tot);



	return out;

}

double  sc::recombM(vector< vector<eccs> > fix2, pars P1) {



	double vnvm=0;
	double totm4=0;
	for (int i=0;i<10;i++){
		SCsub fix=SCM(fix2[i], P1);
		vnvm+=fix.vnvm*fix.m4;
		totm4+=fix.m4;
		fixsc[i]=fix.vnvm*fix.m4;
		fixM4[i]=fix.m4;

	}


	vnvm/=totm4;

	//cout << "tots " << totm4 << endl;
	return vnvm;

}







double  sc::recombMfix(vector< vector<eccs> > fix2, pars P1,int j) {



	double vnvm=0;
	double totm4=0;
	for (int i=0;i<10;i++){
		if (i==j) {
		SCsub fix=SCM(fix2[i], P1);
		vnvm+=fix.vnvm*fix.m4;
		totm4+=fix.m4;}
		else {
		vnvm+=fixsc[i];
		totm4+=fixM4[i];
		}
	}


	vnvm/=totm4;
	return vnvm;

}




void sc::runEP(string name, pars p1, vector< vector<eccs> > sall,int j,  std::vector<cent> cens){



	vector<serr> sc32all;
	int csub=0;
	//for (int j=10;j<=60;j=j+10){

	vector< vector<eccs> > sub(&sall[j-10],&sall[j]);
	sc32all.push_back(sc::jackEP (sub,p1) );
	sc32all[csub].cen=j-5;
	//cout << sc32all[csub].cen << " " << sc32all[csub].ans << " " << sc32all[csub].err << endl;
	csub++;
	//}

	int first=1;
	if (j>10) first=0;
	print(name,sc32all,j,first,cens);

}

void sc::runvn(string name, pars p1, vector< vector<eccs> > sall,int j){



	vector<serr> sc32all;
	int csub=0;

	vector< vector<eccs> > sub(&sall[j-5],&sall[j]);
	sc32all.push_back(sc::jackvn (sub,p1) );
	sc32all[csub].cen=j-2.5;
	csub++;


	int first=1;
	if (j>5) first=0;
	print2(name,sc32all,j,first);

}

serr sc::jackEP (  vector< vector<eccs> > sall, pars P1 ){


	serr out;
	out.ans=sc::recombMEP(sall,P1);

	double m=0,tot=0;
	for(int j=0;j<10;j++){
	int len=sall[j].size();
	tot+=len;
	for(int i=0;i<len;i++){
	  vector< vector<eccs> > sub=sall;
	  sub[j].erase( sub[j].begin() +i );
	  double fsub=sc::recombMfixEP(sub,P1,j);
	  m+=(fsub-out.ans)*(fsub-out.ans);

	}}
	out.err=sqrt((tot-1)*m/tot);

	return out;

}

double  sc::recombMEP(vector< vector<eccs> > fix2, pars P1) {


	SCsub (*bob)(vector<eccs>, pars);

	if (P1.n==2&&P1.m==4&&P1.fn==2&&P1.fm==1) bob=EP224;
	else if (P1.n==2&&P1.m==3&&P1.fn==3&&P1.fm==2) bob=EP22233;
	else if (P1.n==2&&P1.m==4&&P1.fn==4&&P1.fm==2) bob=EP222244;
	else if (P1.n==2&&P1.m==4&&P1.fn==6&&P1.fm==3) bob=EP222222444;
	else if (P1.n==2&&P1.m==3&&P1.p==5&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP235;
	else if (P1.n==2&&P1.m==3&&P1.p==5&&P1.fn==4&&P1.fm==1&&P1.fm==1) bob=EP222235;
	else if (P1.n==1&&P1.m==1&&P1.p==2&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP112;
	else if (P1.n==2&&P1.m==4&&P1.p==6&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP246;
	else if (P1.n==1&&P1.m==2&&P1.p==3&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP123;
	else if (P1.n==3&&P1.m==3&&P1.p==6&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP336;
	else if (P1.n==2&&P1.m==2&&P1.p==6&&P1.fn==2&&P1.fm==1&&P1.fm==1) bob=EP2226;


	double vnvm=0;
	double totm4=0;
	for (int i=0;i<10;i++){
		SCsub fix=(*bob)(fix2[i], P1);
		vnvm+=fix.vnvm*fix.m4;
		totm4+=fix.m4;
		fixsc[i]=fix.vnvm*fix.m4;
		fixM4[i]=fix.m4;
	}


	vnvm/=totm4;
	return vnvm;

}



double  sc::recombMfixEP(vector< vector<eccs> > fix2, pars P1,int j) {

	SCsub (*bob)(vector<eccs>, pars);

	if (P1.n==2&&P1.m==4&&P1.fn==2&&P1.fm==1) bob=EP224;
	else if (P1.n==2&&P1.m==3&&P1.fn==3&&P1.fm==2) bob=EP22233;
	else if (P1.n==2&&P1.m==4&&P1.fn==4&&P1.fm==2) bob=EP222244;
	else if (P1.n==2&&P1.m==4&&P1.fn==6&&P1.fm==3) bob=EP222222444;
	else if (P1.n==2&&P1.m==3&&P1.p==5&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP235;
	else if (P1.n==2&&P1.m==3&&P1.p==5&&P1.fn==4&&P1.fm==1&&P1.fm==1) bob=EP222235;
	else if (P1.n==1&&P1.m==1&&P1.p==2&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP112;
	else if (P1.n==2&&P1.m==4&&P1.p==6&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP246;
	else if (P1.n==1&&P1.m==2&&P1.p==3&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP123;
	else if (P1.n==3&&P1.m==3&&P1.p==6&&P1.fn==1&&P1.fm==1&&P1.fm==1) bob=EP336;
	else if (P1.n==2&&P1.m==2&&P1.p==6&&P1.fn==2&&P1.fm==1&&P1.fm==1) bob=EP2226;

	double vnvm=0;
	double totm4=0;
	for (int i=0;i<10;i++){
		if (i==j) {
		SCsub fix=(*bob)(fix2[i], P1);
		vnvm+=fix.vnvm*fix.m4;
		totm4+=fix.m4;}
		else {
		vnvm+=fixsc[i];
		totm4+=fixM4[i];
		}
	}


	vnvm/=totm4;
	return vnvm;

}


serr sc::jackvn (  vector< vector<eccs> > sall, pars P1 ){



	serr out=sc::recombMvn(sall,P1);
	double m1=0,m2=0,tot=0;
	for(int j=0;j<5;j++){
	int len=sall[j].size();
	tot+=len;
	for(int i=0;i<len;i++){
	  vector< vector<eccs> > sub=sall;
	  sub[j].erase( sub[j].begin() +i );
	  serr fsub=sc::recombMfixvn(sub,P1,j);
	  m1+=(fsub.ans-out.ans)*(fsub.ans-out.ans);
	  m2+=(fsub.ans2-out.ans2)*(fsub.ans2-out.ans2);
	}}
	out.err=sqrt((tot-1)*m1/tot);
	out.err2=sqrt((tot-1)*m2/tot);

	return out;

}

serr  sc::recombMvn(vector< vector<eccs> > fix2, pars P1) {


	SCsub (*bob)(vector<eccs>, pars);

	if (P1.n==2&&P1.m==100) bob=v6;


	double v6o=0,v4o=0;
	double totm6=0,totm4=0;
	for (int i=0;i<5;i++){
		SCsub fix=(*bob)(fix2[i], P1);
		v6o+=fix.v6*fix.m6;
		v4o+=fix.v4*fix.m4;
		totm4+=fix.m4;
		totm6+=fix.m6;
		fixv6[i]=fix.v6*fix.m6;
		fixv4[i]=fix.v4*fix.m4;
		fixM4[i]=fix.m4;
		fixM6[i]=fix.m6;
	}

	serr out;
	out.ans=pow(-v4o/totm4,0.25);
	out.ans2=pow(0.25*v6o/totm6,1./6.);
	return out;

}



serr  sc::recombMfixvn(vector< vector<eccs> > fix2, pars P1,int j) {

	SCsub (*bob)(vector<eccs>, pars);

	if (P1.n==2&&P1.m==100) bob=v6;

	double v6o=0,v4o=0;
	double totm6=0,totm4=0;
	for (int i=0;i<5;i++){
		if (i==j) {
		SCsub fix=(*bob)(fix2[i], P1);
		v6o+=fix.v6*fix.m6;
		v4o+=fix.v4*fix.m4;
		totm4+=fix.m4;
		totm6+=fix.m6;
		}
		else {
		totm4+=fixM4[i];
		totm6+=fixM6[i];

		v6o+=fixv6[i];
		v4o+=fixv4[i];
		}
	}

	serr out;
	out.ans=pow(-v4o/totm4,0.25);
	out.ans2=pow(0.25*v6o/totm6,1./6.);
	return out;

}










SCsub sc::SCM(vector<eccs> vec, pars p1){

	SCsub out;
	if (p1.n!=p1.m){
	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,m2all=0,m4all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m4=m2*(vec[i].M-2.)*(vec[i].M-3.);
	double svn2=vec[i].v[p1.n]*vec[i].v[p1.n];
	double svm2=vec[i].v[p1.m]*vec[i].v[p1.m];

	vn2+=svn2*m2;
	vm2+=svm2*m2;
	vnvm+=svn2*svm2*m4;
	m2all+=m2;
	m4all+=m4;

	//cout << vec[i].v[p1.n] << " " << vec[i].v[p1.m] << " " << vec[i].M << endl;
	}



	out.vnvm=vnvm/m4all/(vn2/m2all*vm2/m2all)-1;
	out.m4=m4all;

//		if(p1.n==4&&p1.m==2&& out.vnvm>1){ cout << "SCM >1 "<< p1.n << " " << p1.m << " " <<    out.vnvm << endl;
////		for (int i=0;i<len;i++){
////		cout <<vec[i].M << "  " << vec[i].v[2] << " " << vec[i].v[3] << " " << vec[i].v[4] << endl;
////		}
//		//getchar();
//		}

	}
	else out=v4cum(vec,p1);

	return out;

}

SCsub sc::v4cum(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vn2=0,vn4=0,m2all=0,m4all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m4=m2*(vec[i].M-2.)*(vec[i].M-3.);

	vn2+=pow(vec[i].v[p1.n],2)*m2;
	vn4+=pow(vec[i].v[p1.n],4)*m4;
	m2all+=m2;
	m4all+=m4;
	}

	SCsub out;

	out.vnvm=(2.*vn2/m2all*vn2/m2all-vn4/m4all)/pow(vn2/m2all,2);
	out.m4=m4all;

	return out;

}


//note this is ONLY for NEXUS, must change for mckln with p1.fn*p1.n within the cosine
SCsub sc::EP3M3(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,vp2=0,m2all=0,m3all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	vn2+=pow(vec[i].v[p1.n],2*p1.fn)*m2;
	vm2+=pow(vec[i].v[p1.m],2*p1.fm)*m2;
	vp2+=pow(vec[i].v[p1.p],2*p1.fp)*m2;
	vnvm+=pow(vec[i].v[p1.n],p1.fn)*pow(vec[i].v[p1.m],p1.fm)*pow(vec[i].v[p1.p],p1.fp)*cos(p1.fn*vec[i].psi[p1.n]+p1.fm*vec[i].psi[p1.m]-p1.fp*vec[i].psi[p1.p])*m3;
	m2all+=m2;
	m3all+=m3;
	}


	//double out=vnvm/m3all;

	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m3all/sqrt(vn2*vm2*vp2/pow(m2all,3));
	else out.vnvm=vnvm/m3all;
	out.m4=m3all;

	return out;

}


double sc::v2M(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vn2=0,m2all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	vn2+=vec[i].v[p1.n]*vec[i].v[p1.n]*m2;
	m2all+=m2;
	}

	double out=sqrt(vn2/m2all);


	return out;

}

double sc::v2(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vn2=0;
	for (int i=0;i<len;i++){
	vn2+=vec[i].v[p1.n]*vec[i].v[p1.n];
	}

	double out=sqrt(vn2/len);


	return out;

}

double sc::v23M(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vn2=0,vm2=0,m2all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	vn2+=vec[i].v[p1.n]*vec[i].v[p1.n]*m2;
	vm2+=vec[i].v[p1.m]*vec[i].v[p1.m]*m2;
	m2all+=m2;
	}

	double out=sqrt(vn2/m2all)*sqrt(vm2/m2all);


	return out;

}

double sc::v23(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vn2=0,vm2=0;
	for (int i=0;i<len;i++){
	vn2+=vec[i].v[p1.n]*vec[i].v[p1.n];
	vm2+=vec[i].v[p1.m]*vec[i].v[p1.m];
	}


	double out=sqrt(vn2/len)*sqrt(vm2/len);


	return out;

}




#endif
