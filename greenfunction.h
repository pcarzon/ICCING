#ifndef GreensFunction_H
#define GreensFunction_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <omp.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

namespace GreensFunctions{

    // EVALUATION MACRO FOR GSL INTERPOLATION FUNCTIONS //
    #define EVALUATE_GSL_INTERPOLATOR_2D(Interpolator,xValue,yValue,xAccelerator,yAccelerator,xMinValue,xMaxValue,yMinValue,yMaxValue)  \
    if((xValue)<(xMinValue) || (xValue)>(xMaxValue) || (yValue)<(yMinValue) || (yValue)>(yMaxValue)){\
        if((yValue)>=(yMinValue) && (yValue)<=(yMaxValue)){ \
            if((xValue)<(xMinValue)){return gsl_spline2d_eval(Interpolator,xMinValue,yValue,xAccelerator,yAccelerator);} \
            else{return 0.0;} \
        }\
        else if((yValue)>(yMaxValue)){\
            if((xValue)<(xMinValue)){return gsl_spline2d_eval(Interpolator,xMinValue,yMaxValue,xAccelerator,yAccelerator);} \
	    else{return gsl_spline2d_eval(Interpolator,xValue,yMaxValue,xAccelerator,yAccelerator);} \
        } \
        else{\
	    std::cerr << "#WARNING " << xValue << " " << xMinValue  << " " << xMaxValue << " " << yValue << " " << yMinValue  << " " << yMaxValue  << std::endl; \
            return 0.0;\
        }\
    }\
    else{\
	    return gsl_spline2d_eval(Interpolator,xValue,yValue,xAccelerator,yAccelerator);\
    }\

    ///////////////////////////////////////////////
    // EVALUTE GREENS FUNCTIONS BY INTERPOLATION //
    ///////////////////////////////////////////////

    // EVALUATE BETWEEN dXdTMin,dXdTMax AND wTMin,wTMax //
    extern double dXdTMin,dXdTMax,wTMin,wTMax;

    // GSL INTERPOLATION OBJECTS //
    extern gsl_interp_accel **FsswTAcc,**GsswTAcc;
    extern gsl_interp_accel **FssdXdTAcc,**GssdXdTAcc;

    extern gsl_spline2d *FssInt,*GssInt;

    extern double FssScalingCurve(double wT,double dXdT);

    extern double Fss(double wT,double dXdT);


    extern double GssScalingCurve(double wT,double dXdT);

    extern double Gss(double wT,double dXdT);


    /////////////////////////////////////////////////////////
    //    INTERPOLATE COORDINATE SPACE GRRENS FUNCTIONS    //
    // AS FUNCTIONS OF wTilde and (\Delta x)/(\Delta \tau) //
    /////////////////////////////////////////////////////////

    extern double *wTValues,*dXdTValues;

    extern double *FssValues;
    extern double *GssValues;

    extern void Setup(int NumberOfTimes,int NumberOfPoints);
    // READ INPUT FILE //
    // INPUT FILE MUST HAVE FOLLOWING STRUCTURE: 1:wTilde 2:|x-x_0|/|tau-tau_0| 3:|tau-tau_0|^2*Fss 4:|tau-tau_0|^2*Gss //
    // AVOID EMPTY LINES IN INPUT FILE! //
    extern void SetValues(std::string fname,int NumberOfTimes,int NumberOfPoints);

    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////

    extern void SetupInterpolators(int NumberOfTimes,int NumberOfPoints);


    // CREATE OUTPUT //
    extern void Output(std::string fname,int NwT,int NdXdT);


} // GreensFunctions
#endif
