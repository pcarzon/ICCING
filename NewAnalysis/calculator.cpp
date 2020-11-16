#include "calculator.h"

Calculator::Calculator(vector<vector<Event>> sorted_events)
{
  binned_events = sorted_events;
  for (int bin = 0; bin < binned_events.size(); bin++)
  {
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //	Calculate Values
    int evs = binned_events[bin].size();
  /*  for(int ev = 0; ev < evs; ev++)
    {
      //	Calculate average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        avg_Vn[n] += pow(binned_events[bin][ev].eccentricity[n], 2)/evs;
        avg_Vn_2nd[n] += pow(binned_events[bin][ev].eccentricity[n], 4)/evs;
        avg_Vn_3rd[n] += pow(binned_events[bin][ev].eccentricity[n], 6)/evs;
      }
    }
*/
  }
}
//__________________________________________________________________________________________
//##########################################################################################
//  Class Destructor
//##########################################################################################
Calculator::~Calculator()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Calculator Copy Function
//##########################################################################################
void Calculator::CopyCalculator(const Calculator &e)
{
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Calculator::Calculator(const Calculator &original)
{
  CopyCalculator(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Calculator& Calculator::operator= (const Calculator& original)
{
	CopyCalculator(original);
	return *this;
}
//__________________________________________________________________________________________

void Calculator::CalculateCummulants()
{
	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate Values
		int evs = binned_events[bin].size();
		double Vn_4part[7] = {0}, Vn_6part[7] = {0};
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
        EvErr_avg_Vn[n] = (avg_Vn[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
        EvErr_avg_Vn_2nd[n] = (avg_Vn_2nd[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 4))/(evs - 1);
        EvErr_avg_Vn_3rd[n] = (avg_Vn_3rd[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 6))/(evs - 1);
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

	}

}

void Calculator::Calculate_2Particle_Cummulants()
{
  cout << "2 particle cummulants " << binned_events.size() << endl;

	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();
		double Vn_2part[7] = {0};
//cout << "1" << endl;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values

    for(int ev = 0; ev < evs; ev++)
    {
      //	Calculate average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        avg_Vn[n] += pow(binned_events[bin][ev].eccentricity[n], 2)/evs;
        avg_Vn_2nd[n] += pow(binned_events[bin][ev].eccentricity[n], 4)/evs;
        avg_Vn_3rd[n] += pow(binned_events[bin][ev].eccentricity[n], 6)/evs;
      }
    }

		//	Calculate 2, 4, and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
  //    cout << "2" << endl;
			Vn_2part[n] = sqrt(avg_Vn[n]);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
    double EvErr_avg_Vn[7] = {0};
		double err_Vn_2part[7] = {0};

		for(int ev = 0; ev < evs; ev++)
		{
//      cout << "3" << endl;
      // Calculate error by event of average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        EvErr_avg_Vn[n] = (avg_Vn[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
      }

			//	Calculate error for 2, 4, and 6 particle cummulants
		 	for (int n = 2; n <= 5; n++)
			{
		 		err_Vn_2part[n] += pow(Vn_2part[n] - sqrt(EvErr_avg_Vn[n]), 2);
			}
		}

		//	Normalize errors
		for (int n = 2; n <= 5; n++)
		{
			err_Vn_2part[n] = sqrt(err_Vn_2part[n]*(evs - 1)/evs);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //  cout << "Got here" << endl;
    cout << endl << (0.5 + bin)*(100/binned_events.size()) << " ";

    for (int n = 2; n <= 5; n++)
    cout <<  Vn_2part[n] << " " << err_Vn_2part[n] << " ";

	}

}

void Calculator::Calculate_4and6Particle_Cummulants()
{
	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();
		double Vn_4part[7] = {0}, Vn_6part[7] = {0};

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values

		//	Calculate 2, 4, and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
			Vn_4part[n] = 2.*pow(avg_Vn[n], 2) - avg_Vn_2nd[n];
			Vn_6part[n] = 0.25*(avg_Vn_3rd[n] - 9.*avg_Vn[n]*avg_Vn_2nd[n] + 12.*pow(avg_Vn[n], 3));
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
    double EvErr_avg_Vn[7] = {0}, EvErr_avg_Vn_2nd[7] = {0},  EvErr_avg_Vn_3rd[7] = {0};
	  double err_Vn_4part[7] = {0},  err_Vn_6part[7] = {0};

		for(int ev = 0; ev < evs; ev++)
		{
      // Calculate error by event of average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        EvErr_avg_Vn[n] = (avg_Vn[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
        EvErr_avg_Vn_2nd[n] = (avg_Vn_2nd[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 4))/(evs - 1);
        EvErr_avg_Vn_3rd[n] = (avg_Vn_3rd[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 6))/(evs - 1);
      }

			//	Calculate error for 2, 4, and 6 particle cummulants
		 	for (int n = 2; n <= 5; n++)
			{
				err_Vn_4part[n] += pow(Vn_4part[n] - (2.*pow(EvErr_avg_Vn[n], 2) - EvErr_avg_Vn_2nd[n]), 2);
				err_Vn_6part[n] += pow(Vn_6part[n] - 0.25*(EvErr_avg_Vn_3rd[n] - 9.*EvErr_avg_Vn[n]*EvErr_avg_Vn_2nd[n] + 12.*pow(EvErr_avg_Vn[n], 3)), 2);
			}

		}

		//	Normalize errors
		for (int n = 2; n <= 5; n++)
		{
			err_Vn_4part[n] = sqrt(err_Vn_4part[n]*(evs - 1)/evs);
			err_Vn_6part[n] = sqrt(err_Vn_6part[n]*(evs - 1)/evs);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	}

}

void Calculator::Calculate_V2_V3_CummulantRatio()
{
	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();
		double V2_2part_V3_2part = 0;

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values

		V2_2part_V3_2part = sqrt(avg_Vn[2]/avg_Vn[3]);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
		double EvErr_avg_Vn[7] = {0};
    double err_V2_2part_V3_2part = 0;

		for(int ev = 0; ev < evs; ev++)
		{
      // Calculate error by event of average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        EvErr_avg_Vn[n] = (avg_Vn[n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
      }

			err_V2_2part_V3_2part += pow(V2_2part_V3_2part - sqrt(EvErr_avg_Vn[2]/EvErr_avg_Vn[3]), 2);
		}

		//	Normalize errors
		err_V2_2part_V3_2part = sqrt(err_V2_2part_V3_2part*(evs - 1)/evs);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	}

}

void Calculator::CalculatePtObservables()
{
	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();

		double mean_pt = 0, mean_pt_2nd = 0, RMS_pt = 0;

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values
		for(int ev = 0; ev < evs; ev++) mean_pt += binned_events[bin][ev].average_pt;
		mean_pt /= evs;

      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //	Calculate Values
      for(int ev = 0; ev < evs; ev++)
      {
        //	Calculate average Vn, Vn^2, Vn^3
        mean_pt_2nd += pow(binned_events[bin][ev].average_pt - mean_pt, 2.);
      }


		mean_pt_2nd /= evs;
		RMS_pt = sqrt(mean_pt_2nd)/mean_pt;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
		double err_mean_pt = 0, err_mean_pt_2nd = 0, err_RMS_pt = 0;

		for(int ev = 0; ev < evs; ev++)
		{
			double EvErr_mean_pt = (mean_pt*evs - binned_events[bin][ev].average_pt)/(evs - 1.);
			double EvErr_mean_pt_2nd = (mean_pt_2nd*evs - pow(binned_events[bin][ev].average_pt - EvErr_mean_pt, 2))/(evs - 1);

		 	err_mean_pt += pow(mean_pt - EvErr_mean_pt, 2);
		 	err_mean_pt_2nd += pow(mean_pt_2nd - EvErr_mean_pt_2nd, 2);
		 	err_RMS_pt = pow(RMS_pt - sqrt(EvErr_mean_pt_2nd)/EvErr_mean_pt, 2);
		}

		err_mean_pt = sqrt(err_mean_pt*(evs - 1)/evs);
		err_mean_pt_2nd = sqrt(err_mean_pt_2nd*(evs - 1)/evs);
		err_RMS_pt = sqrt(err_RMS_pt*(evs - 1)/evs);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	}

}
