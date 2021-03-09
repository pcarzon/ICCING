#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "BackgroundAttractor.cpp"
#include "GreensFunctions.cpp"

int main(){
    
    ///////////////////////////////////////////////////////////////////////////////
    // SET VALUES FOR BACKGROUND ATTRACTOR AND INTERPOLATION OF GREENS FUNCTIONS //
    ///////////////////////////////////////////////////////////////////////////////

    // SET NUMBER OF POINTS FOR BACKGROUND ATTRACTOR //
    int NumberOfPointsBackground=386;

    // SET NUMBER OF dX/dTAU-POINTS AND NUMBER OF wTilde-POINTS FOR INTERPOLATION OF GREENS FUNCTIONS //
    int NumberOfPointsGreensFunctions=128;
    int NumberOfTimesGreensFunctions=100;



    // SET VALUE OF EtaOverS, C_\infty AN (e\tau)_0 //
    double EtaOverS=1/(4.0*M_PI);
    double CInfinity=0.800226;
    double eTau0=1.121358171;



    //////////////////////////////////
    // COMPUTE BACKGROUND ATTRACTOR //
    //////////////////////////////////

    // GET VALUES FROM INPUT FILE AND SETUP INTERPOLATION OF BACKGROUND ATTRACTOR //
    BackgroundAttractor::Setup(NumberOfPointsBackground,"BackgroundAttractorRTA.txt",CInfinity);

    // CREATE OUTPUT FOR BACKGROUND ATTRACTOR //
    BackgroundAttractor::Output("BackgroundEnergy.txt",EtaOverS,eTau0);




    ///////////////////////////////////////////////////////////////////////////////
    // COMPUTE INTERPOLATION OF GREENS FUNCTIONS FOR CHARGES IN COORDINATE SPACE //
    ///////////////////////////////////////////////////////////////////////////////

    // SETUP INTERPOLATION //
    GreensFunctions::Setup(NumberOfTimesGreensFunctions,NumberOfPointsGreensFunctions);
    
    // GET VALUE FROM INPUT FILE AND SET VALUES FOR INTERPOLATION //
    GreensFunctions::SetValues("ChargeGreensFunctionsRTA.txt",NumberOfTimesGreensFunctions,NumberOfPointsGreensFunctions);
    
    // INITIALIZE THE INTERPOLATION //
    GreensFunctions::SetupInterpolators(NumberOfTimesGreensFunctions,NumberOfPointsGreensFunctions);
    
    // EVALUATE INTERPOLATION AND CREATE OUTPUT //
    GreensFunctions::Output("InterpolatedChargeGreensFunctions.txt",1000,1000);
    
    

    
}
