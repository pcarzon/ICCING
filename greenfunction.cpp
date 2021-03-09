#include "greenfunction.h"

namespace GreensFunctions{

    ///////////////////////////////////////////////
    // EVALUTE GREENS FUNCTIONS BY INTERPOLATION //
    ///////////////////////////////////////////////
    // EVALUATE BETWEEN dXdTMin,dXdTMax AND wTMin,wTMax //
    double dXdTMin,dXdTMax,wTMin,wTMax;

    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **FsswTAcc,**GsswTAcc;
    gsl_interp_accel **FssdXdTAcc,**GssdXdTAcc;

    gsl_spline2d *FssInt,*GssInt;


    double FssScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(FssInt,wT,dXdT,FsswTAcc[tID],FssdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
    } // FssScalingCurve

    double Fss(double wT,double dXdT){
        return FssScalingCurve(wT,dXdT);
    } // Fss


    double GssScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GssInt,wT,dXdT,FsswTAcc[tID],GssdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
    } // GssScalingCurve

    double Gss(double wT,double dXdT){
        return GssScalingCurve(wT,dXdT);
    } // Gss


    /////////////////////////////////////////////////////////
    //    INTERPOLATE COORDINATE SPACE GRRENS FUNCTIONS    //
    // AS FUNCTIONS OF wTilde and (\Delta x)/(\Delta \tau) //
    /////////////////////////////////////////////////////////
    double *wTValues,*dXdTValues;

    double *FssValues;
    double *GssValues;


    void Setup(int NumberOfTimes,int NumberOfPoints){

        // ALLOCATE //
        wTValues=new double[NumberOfTimes];
        dXdTValues=new double[NumberOfPoints];

        FssValues=new double[NumberOfTimes*NumberOfPoints];
        GssValues=new double[NumberOfTimes*NumberOfPoints];

        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();

        FsswTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GsswTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];

        FssdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GssdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];

        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){

            FsswTAcc[i] = gsl_interp_accel_alloc ();
            GsswTAcc[i] = gsl_interp_accel_alloc ();

            FssdXdTAcc[i] = gsl_interp_accel_alloc ();
            GssdXdTAcc[i] = gsl_interp_accel_alloc ();
            }

        FssInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        GssInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);

    } // Setup

    // READ INPUT FILE //
    // INPUT FILE MUST HAVE FOLLOWING STRUCTURE: 1:wTilde 2:|x-x_0|/|tau-tau_0| 3:|tau-tau_0|^2*Fss 4:|tau-tau_0|^2*Gss //
    // AVOID EMPTY LINES IN INPUT FILE! //
    void SetValues(std::string fname,int NumberOfTimes,int NumberOfPoints){

        double FssVal[NumberOfTimes*NumberOfPoints];
        double GssVal[NumberOfTimes*NumberOfPoints];

        std::ifstream InStream;
        InStream.open(fname);

        int wCounter=0; int xCounter=0;

        // READ INPUT FILE //
        while(InStream.good()){

            double wT; double dXdT; double Fss; double Gss;

            InStream >> wT; InStream >> dXdT; InStream >> Fss; InStream >> Gss;


            // WRITE EACH POSITION-STEP ONCE INTO dXdTValues //
            dXdTValues[xCounter]=dXdT;

            // WRITE EACH TIME-STEP ONCE INTO wTValues //
            if((wCounter%NumberOfPoints)==0){
                wTValues[wCounter/NumberOfPoints]=wT;
            }

            // WRITE VALUES OF GREENS FUNCTIONS INTO ARRAYS //
            FssVal[wCounter]=Fss; GssVal[wCounter]=Gss;


            wCounter++;

            // RESET xCounter WHEN NEW BLOCK IS REACHED IN INPUT FILE //
            if(xCounter<NumberOfPoints){
                xCounter++;
            }

            if(xCounter==NumberOfPoints){
                xCounter=0;
            }

        }

        // CHECK WHETHER TIMES wTilde AND POSITIONS dXdT ARE READ CORRECTLY //
        /* int k=0;
        while (k < NumberOfTimes) {
           std::cout << wTValues[k] << " ";
           printf("\n");
           k++;
        }

        int l=0;
        while (l < NumberOfPoints) {
           std::cout << dXdTValues[l] << " ";
           printf("\n");
           l++;
        } */


        // SET GRID VALUES FOR INTERPOLATION //
        int FIndex=0;
        for(int wTIndex=0;wTIndex<NumberOfTimes;wTIndex++){

            for(int xIndex=0;xIndex<NumberOfPoints;xIndex++){

                gsl_spline2d_set(FssInt,FssValues,wTIndex,xIndex,FssVal[FIndex]);
                gsl_spline2d_set(GssInt,GssValues,wTIndex,xIndex,GssVal[FIndex]);

                FIndex++;

            }

        }

    } // SetValues

    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////

    void SetupInterpolators(int NumberOfTimes,int NumberOfPoints){

        // SET BOUNDARIES //
        wTMin=wTValues[0]; wTMax=wTValues[NumberOfTimes-1];
        dXdTMin=dXdTValues[0]; dXdTMax=dXdTValues[NumberOfPoints-1];

        // INITIALIZE INTERPOLATOR //
        gsl_spline2d_init(FssInt,wTValues,dXdTValues,FssValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(GssInt,wTValues,dXdTValues,GssValues,NumberOfTimes,NumberOfPoints);

        // CLEAN-UP //
        delete[] wTValues;
        delete[] dXdTValues;

        delete[] FssValues;
        delete[] GssValues;

    } // SetupInterpolators


    // CREATE OUTPUT //
    void Output(std::string fname,int NwT,int NdXdT){


        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# 1:wTilde=tau T(tau)/(eta/s) 2:|x-x_0|/|tau-tau0| 3:|tau-tau0|^2Fss 4:|tau-tau0|^2Gss" << std::endl;

        for(int tIndex=0;tIndex<NwT;tIndex++){

            double wTi=wTMin + (tIndex+0.5)*(wTMax-wTMin)/(NwT);

            for(int xIndex=0;xIndex<NdXdT;xIndex++){

                double dXdTi=dXdTMin + (xIndex+0.5)*(dXdTMax-dXdTMin)/(NdXdT);

                Outstream << wTi << " " << dXdTi << " " << Fss(wTi,dXdTi) << " " << Gss(wTi,dXdTi) << std::endl;

            }

            Outstream << std::endl;
            Outstream << std::endl;

        }

        Outstream.close();
    } // Output


} // GreensFunctions
