#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace BackgroundAttractor {

    // SET DEGREES OF FREEDOM //
    double Nc=3.0;
    double Nf=3.0;
    double nuG=2.0*(Nc*Nc-1.0);
    double nuQ=2.0*Nc*Nf; 
    double nuEff=nuG+7.0/4.0*nuQ;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **EAcc;
    gsl_spline *EInt;
    
    double wTMin; double wTMax; double CInfty;
    
    // SET-UP VALUES FOR wTilde (wT) AND ENERGY ATTRACTOR (EVal) FORM INPUT FILE //
    // INPUT FILE MUST HAVE FOLLOWING STRUCTURE: 1:wTilde 2:EnergyAttractor //
    // AVOID EMPTY LINES IN INPUT FILE! //
    void Setup(int NumberOfPoints,std::string fname,double CInftyVal){
        
        // SET C_\infty //
        CInfty = CInftyVal;

        // SET DATA //
        double wTildeValues[NumberOfPoints];
        double EValues[NumberOfPoints];
        
        std::ifstream InStream;
        InStream.open(fname);
        
        int i=0;
        
        while(InStream.good()){
            
            double wT; double EVal;
            
            InStream >> wT; InStream >> EVal;
            
            wTildeValues[i]=wT; EValues[i]=EVal;
            
            i++;
            
        }

        // CHECK WHETHER TIMES wTilde AND ENERGY ATTRACTOR E ARE READ CORRECTLY //
        /* int k=0;
        while (k < NumberOfPoints) {
           std::cout << wTildeValues[k] << " " << EValues[k] << " ";
           printf("\n");
           k++;
        } */

        // SETUP SPLINE //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        EAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];

        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            EAcc[i] = gsl_interp_accel_alloc ();
        }

        EInt=gsl_spline_alloc(gsl_interp_cspline,NumberOfPoints);
        gsl_spline_init(EInt,wTildeValues,EValues,NumberOfPoints);

        // SET BOUNDARIES //
        wTMin=wTildeValues[0];
        wTMax=wTildeValues[NumberOfPoints-1];

        
    } // Setup
    
    // ENERGY ATTRACTOR CURVE //
    double E(double wT){
        
        if(wT<wTMin){
            return 1.0/CInfty*std::pow(wT,4.0/9.0);
        }
        else if(wT>wTMax){
            return 1.0-2.0/(3.0*M_PI*wT);
        }
        else{
            int tID=omp_get_thread_num();
            return gsl_spline_eval(EInt,wT,EAcc[tID]);
        }
        
    } // E

    void GetValues(double eTau0,double Tau,double etaOverS,double &e,double &wTilde){

        // DETERMINE (e(tau) tau^{4/3})_{infty} //
        double eTau43Infty=std::pow(4.0*M_PI*etaOverS,4.0/9.0)*std::pow(M_PI*M_PI*nuEff/30.0,1.0/9.0)*CInfty*std::pow(eTau0,8.0/9.0);       

        //////////////////////////////////////////////////////////
        // DETERMINE TEMPERATURE SELF-CONSISTENTLY ACCORDING TO //
        // e(T)tau^{4/3} = E(wTilde) (e(tau) tau^{4/3})_{infty} //
        //          wTilde= (T tau)/(4pi eta/s)                 //
        //////////////////////////////////////////////////////////
        
        double TLow=0.0; double THigh=std::pow(eTau43Infty/((M_PI*M_PI/30.0)*nuEff*std::pow(1.0,4.0)*std::pow(Tau,4.0/3.0)),1.0/4.0);
        
        double TMid=(THigh+TLow)/2.0;
        double wTildeMid=(TMid*Tau)/(4.0*M_PI*etaOverS);


        while(THigh-TLow>1E-6*TMid){
            
            if(E(wTildeMid)/std::pow(TMid,4)>(M_PI*M_PI/30.0)*nuEff*std::pow(1.0,4.0)*std::pow(Tau,4.0/3.0)/eTau43Infty){
                TLow=TMid;
            }
            else{
                THigh=TMid;
            }
            
            TMid=(THigh+TLow)/2.0;
            wTildeMid=(TMid*Tau)/(4.0*M_PI*etaOverS);
            
        }
        
        // CHECK THAT eEq(T) == E(wTilde) (e(tau) tau^{4/3})_{infty} IS SOLVED //
        //std::cerr << "wT=" << wTilde << " " << "eEq=" << (M_PI*M_PI/30.0)*nuEff*std::pow(TMid,4.0) << " " << "e=" << eTau43Infty*E((TMid*Tau)/(4.0*M_PI*etaOverS))/std::pow(Tau,4.0/3.0) << std::endl;
        
        // SET FINAL VALUE OF wTilde //
        wTilde=(TMid*Tau)/(4.0*M_PI*etaOverS);
        
        // SET FINAL VALUES OF T,e IN GeV //
        e=(M_PI*M_PI/30.0)*nuEff*std::pow(TMid,4.0);



    } // GetValues


    // CREATE OUTPUT //    
    void Output(std::string fname,double EtaOverS,double eTau0){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        
        Outstream << "# 1:Tau 2:wTilde=tau T(tau)/(eta/s) 3:e" << std::endl;
        
        for(int i=0;i<4096;i++){
        
            double Tau=0.0001*exp(0.005*i);
        
            // DETERMINE e(tau) and wTilde(tau) BASED ON BACKGROUND-ATTRACTOR //
            double e,wTilde;
            GetValues(eTau0,Tau,EtaOverS,e,wTilde);

            Outstream << Tau << " " << wTilde << " " << e << std::endl;
        }

        Outstream << std::endl;

        Outstream.close();
    } // Output
    
} // BackgroundAttractor

