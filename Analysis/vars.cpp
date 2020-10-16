#ifndef _VARS_CPP_
#define _VARS_CPP_

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

# include "vars.h"


using namespace std;



SCsub EP224(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vnn2=0,vm2=0,m2all=0,m3all=0,m4all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	vn2+=pow(vec[i].v[2],2)*m2;
	vm2+=pow(vec[i].v[4],2)*m2;
	vnn2+=pow(vec[i].v[2],4)*m4;
	vnvm+=pow(vec[i].v[2],2)*vec[i].v[4]*cos(4.*vec[i].psi[2]-4.*vec[i].psi[4])*m3;
	//vnvm+=pow(vec[i].v[2],2)*vec[i].v[4]*cos(4*vec[i].psi[2]-4*vec[i].psi[4])*m3;
	m2all+=m2;
	m3all+=m3;
	m4all+=m4;
	}
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m3all/sqrt(vn2*vn2*vm2/pow(m2all,3));
	else if (p1.swit==2) out.vnvm=vnvm/m3all/sqrt(vnn2*vm2/(m2all*m4all));
	else out.vnvm=vnvm/m3all;
	out.m4=m3all;
	
	return out;

}

SCsub EP2226(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vnn2=0,vm2=0,m2all=0,m4all=0,m6all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	double m6=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.);
	vn2+=pow(vec[i].v[2],2)*m2;
	vm2+=pow(vec[i].v[6],2)*m2;
	vnn2+=pow(vec[i].v[2],6)*m6;
	vnvm+=pow(vec[i].v[2],3)*vec[i].v[6]*cos(6.*vec[i].psi[2]-6.*vec[i].psi[6])*m4;
	//vnvm+=pow(vec[i].v[2],2)*vec[i].v[4]*cos(4*vec[i].psi[2]-4*vec[i].psi[4])*m3;
	m2all+=m2;
	m4all+=m4;
	m6all+=m6;
	}
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m4all/sqrt(vn2*vn2*vn2*vm2/pow(m2all,4));
	else if (p1.swit==2) out.vnvm=vnvm/m4all/sqrt(vnn2*vm2/(m2all*m6all));
	else out.vnvm=vnvm/m4all;
	out.m4=m4all;
	
	//if (p1.swit==2) cout << out.vnvm << " " << m2all << " " << m6all << endl;
	
	return out;

}

SCsub EP22233(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,m4all=0,m5all=0,m6all=0;
	for (int i=0;i<len;i++){
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	double m5=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.);
	double m6=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.);
	vn2+=pow(vec[i].v[2],2*3)*m6;
	vm2+=pow(vec[i].v[3],2*2)*m4;
	vnvm+=pow(vec[i].v[2],3)*pow(vec[i].v[3],2)*cos(3*vec[i].psi[2]-2*vec[i].psi[3])*m5;
	m4all+=m4;
	m5all+=m5;
	m6all+=m6;
	}
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m5all/sqrt(vn2*vm2/(m4all*m6all));
	else out.vnvm=vnvm/m5all;
	out.m4=m5all;
	
	return out;

}

SCsub EP222244(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,m4all=0,m6all=0,m8all=0;
	for (int i=0;i<len;i++){
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	double m6=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.);
	double m8=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.)*(vec[i].M-6.)*(vec[i].M-7.);
	vn2+=pow(vec[i].v[2],2*4)*m8;
	vm2+=pow(vec[i].v[4],2*2)*m4;
	vnvm+=pow(vec[i].v[2],4)*pow(vec[i].v[4],2)*cos(4*vec[i].psi[2]-2*vec[i].psi[4])*m6;
	m4all+=m4;
	m6all+=m6;
	m8all+=m8;
	}
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m6all/sqrt(vn2*vm2/(m4all*m8all));
	else out.vnvm=vnvm/m6all;
	out.m4=m6all;
	
	return out;

}

SCsub EP222222444(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,m6all=0,m9all=0,m12all=0;
	for (int i=0;i<len;i++){
	double m6=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.);
	double m9=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.)*(vec[i].M-6.)*(vec[i].M-7.)*(vec[i].M-8.);
	double m12=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.)*(vec[i].M-5.)*(vec[i].M-6.)*(vec[i].M-7.)*(vec[i].M-8.)*(vec[i].M-9.)*(vec[i].M-10.)*(vec[i].M-11.);
	vn2+=pow(vec[i].v[2],2*6)*m12;
	vm2+=pow(vec[i].v[4],2*3)*m6;
	vnvm+=pow(vec[i].v[2],6)*pow(vec[i].v[4],3)*cos(6*vec[i].psi[2]-3*vec[i].psi[4])*m9;
	
	m6all+=m6;
	m9all+=m9;
	m12all+=m12;
	}
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m9all/sqrt(vn2*vm2/(m6all*m12all));
	else out.vnvm=vnvm/m9all;
	out.m4=m6all;
	
	return out;

}

SCsub EP235(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,vp2=0,vnm2=0,m2all=0,m3all=0,m4all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	vn2+=pow(vec[i].v[2],2)*m2;
	vm2+=pow(vec[i].v[3],2)*m2;
	vp2+=pow(vec[i].v[5],2)*m2;
	vnm2+=pow(vec[i].v[2],2)*pow(vec[i].v[3],2)*m4;
	vnvm+=vec[i].v[2]*vec[i].v[3]*vec[i].v[5]*cos(2.*vec[i].psi[2]+3.*vec[i].psi[3]-5.*vec[i].psi[5])*m3;
	m2all+=m2;
	m3all+=m3;
	m4all+=m4;
	}
	
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m3all/sqrt(vn2*vm2*vp2/pow(m2all,3));
	else if (p1.swit==2) out.vnvm=vnvm/m3all/sqrt(vnm2*vp2/(m2all*m4all));
	else out.vnvm=vnvm/m3all;
	out.m4=m3all;
	
	return out;

}

SCsub EP112(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,vp2=0,m2all=0,m3all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	vn2+=pow(vec[i].v[p1.n],2)*m2;
	vm2+=pow(vec[i].v[p1.m],2)*m2;
	vp2+=pow(vec[i].v[p1.p],2)*m2;
	vnvm+=pow(vec[i].v[1],2)*vec[i].v[2]*cos(2.*vec[i].psi[1]-2.*vec[i].psi[2])*m3;
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

SCsub EP123(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,vp2=0,m2all=0,m3all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	vn2+=pow(vec[i].v[p1.n],2)*m2;
	vm2+=pow(vec[i].v[p1.m],2)*m2;
	vp2+=pow(vec[i].v[p1.p],2)*m2;
	vnvm+=vec[i].v[1]*vec[i].v[2]*vec[i].v[3]*cos(vec[i].psi[1]+2.*vec[i].psi[2]-3.*vec[i].psi[3])*m3;
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

SCsub EP336(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,vp2=0,m2all=0,m3all=0,m4all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	vn2+=pow(vec[i].v[3],2)*m2;
	vm2+=pow(vec[i].v[6],2)*m2;
	vp2+=pow(vec[i].v[3],4)*m4;
	vnvm+=vec[i].v[p1.n]*vec[i].v[p1.m]*vec[i].v[p1.p]*cos(3*vec[i].psi[p1.n]+3*vec[i].psi[p1.m]-6.*vec[i].psi[p1.p])*m3;
	m2all+=m2;
	m3all+=m3;
	m4all+=m4;
	}
	
	
	//double out=vnvm/m3all;
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m3all/sqrt(vn2*vn2*vm2/pow(m2all,3));
	else if (p1.swit==2) out.vnvm=vnvm/m3all/sqrt(vp2*vm2/(m2all*m4all));
	else out.vnvm=vnvm/m3all;
	out.m4=m3all;
	
	return out;

}

SCsub EP246(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vnm2=0,vn2=0,vm2=0,vp2=0,m2all=0,m3all=0,m4all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.);
	vnm2+=pow(vec[i].v[p1.n]*vec[i].v[p1.m],2)*m4;
	vn2+=pow(vec[i].v[p1.n],2)*m2;
	vm2+=pow(vec[i].v[p1.m],2)*m2;
	vp2+=pow(vec[i].v[p1.p],2)*m2;
	vnvm+=vec[i].v[p1.n]*vec[i].v[p1.m]*vec[i].v[p1.p]*cos(2*vec[i].psi[p1.n]+4*vec[i].psi[p1.m]-6.*vec[i].psi[p1.p])*m3;
	m2all+=m2;
	m3all+=m3;
	m4all+=m4;
	}
	
	
	//double out=vnvm/m3all;
	
	SCsub out;
	if (p1.swit==1) out.vnvm=vnvm/m3all/sqrt(vn2*vm2*vp2/pow(m2all,3));
	else if (p1.swit==2) out.vnvm=vnvm/m3all/sqrt(vp2*vnm2/(m2all*m4all));
	else out.vnvm=vnvm/m3all;
	out.m4=m3all;
	
	return out;

}


SCsub EP222235(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vnvm=0,vn2=0,vm2=0,vp2=0,m2all=0,m3all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m3=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	vn2+=pow(vec[i].v[p1.n],2)*m2;
	vm2+=pow(vec[i].v[p1.m],2)*m2;
	vp2+=pow(vec[i].v[p1.p],2)*m2;
	vnvm+=vec[i].v[2]*vec[i].v[3]*vec[i].v[5]*cos(vec[i].psi[3]+vec[i].psi[5]-4*vec[i].psi[2])*m3;
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

SCsub v6(vector<eccs> vec, pars p1){

	int len=vec.size();
	double vn2=0,vn4=0,vn6=0,m2all=0,m4all=0,m6all=0;
	for (int i=0;i<len;i++){
	double m2=vec[i].M*(vec[i].M-1.);
	double m4=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.);
	double m6=vec[i].M*(vec[i].M-1.)*(vec[i].M-2.)*(vec[i].M-3.)*(vec[i].M-4.);
	vn6+=pow(vec[i].v[p1.n],6)*m6;
	vn4+=pow(vec[i].v[p1.n],4)*m4;
	vn2+=pow(vec[i].v[p1.n],2)*m2;
	m2all+=m2;
	m4all+=m4;
	m6all+=m6;
	}
	
	
	
	double var=vn2/m2all;
	double kur=vn4/m4all;
	
	SCsub out;
	out.v4=kur-2.*pow(var,2);
	out.v6=vn6/m6all-9.*kur*var+12.*pow(var,3);
	out.m4=m4all;
	out.m6=m6all;
	
	return out;

}


#endif


