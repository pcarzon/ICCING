#include "calculator.h"

Calculator::Calculator(vector<vector<Event>> sorted_events)
{
  binned_events = sorted_events;
  avg_Vn.resize(binned_events.size(), vector<double>(7, 0));
  avg_Vn_2nd.resize(binned_events.size(), vector<double>(7, 0));
  avg_Vn_3rd.resize(binned_events.size(), vector<double>(7, 0));

  for (int bin = 0; bin < binned_events.size(); bin++)
  {
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //	Calculate Values
    int evs = binned_events[bin].size();
    for(int ev = 0; ev < evs; ev++)
    {
      //	Calculate average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        avg_Vn[bin][n] += pow(binned_events[bin][ev].eccentricity[n], 2)/evs;
        avg_Vn_2nd[bin][n] += pow(binned_events[bin][ev].eccentricity[n], 4)/evs;
        avg_Vn_3rd[bin][n] += pow(binned_events[bin][ev].eccentricity[n], 6)/evs;
      }
    }

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

vector<vector<double>> Calculator::CalculateCummulants()
{
  vector<vector<double>> Vn_4part, Vn_6part, err_Vn_4part, err_Vn_4part_4th, err_Vn_6part, output;
  Vn_4part.resize(binned_events.size(), vector<double>(7, 0));
  Vn_6part.resize(binned_events.size(), vector<double>(7, 0));
  err_Vn_4part.resize(binned_events.size(), vector<double>(7, 0));
  err_Vn_4part_4th.resize(binned_events.size(), vector<double>(7, 0));
  err_Vn_6part.resize(binned_events.size(), vector<double>(7, 0));

	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate Values
		int evs = binned_events[bin].size();
		//	Calculate 4 and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
			Vn_4part[bin][n] = pow(2.*pow(avg_Vn[bin][n], 2)
												- avg_Vn_2nd[bin][n], 0.25)
										/sqrt(avg_Vn[bin][n]);
			Vn_6part[bin][n] = pow(0.25*(avg_Vn_3rd[bin][n]
													- 9.*avg_Vn[bin][n]*avg_Vn_2nd[bin][n]
													+ 12.*pow(avg_Vn[bin][n], 3)), 1./6.)
										/pow(2.*pow(avg_Vn[bin][n], 2) - avg_Vn_2nd[bin][n], 0.25);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate Errors
    double EvErr_avg_Vn[7] = {0}, EvErr_avg_Vn_2nd[7] = {0},  EvErr_avg_Vn_3rd[7] = {0};
		for(int ev = 0; ev < evs; ev++)
		{

      // Calculate error by event of average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        EvErr_avg_Vn[n] = (avg_Vn[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
        EvErr_avg_Vn_2nd[n] = (avg_Vn_2nd[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 4))/(evs - 1);
        EvErr_avg_Vn_3rd[n] = (avg_Vn_3rd[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 6))/(evs - 1);
      }

			//	Calculate error for 4 and 6 particle cummulants (Also 4 particle cummulant to the 4th power)
		 	for (int n = 2; n <= 5; n++)
			{
		 		err_Vn_4part[bin][n] += pow(Vn_4part[bin][n]
																- pow(2.*pow(EvErr_avg_Vn[n], 2)
																 			- EvErr_avg_Vn_2nd[n], 0.25)
														/sqrt(EvErr_avg_Vn[n]), 2);
      	err_Vn_4part_4th[bin][n] += pow(pow(Vn_4part[bin][n], 4)
																		- (2.*pow(EvErr_avg_Vn[n], 2)
																				- EvErr_avg_Vn_2nd[n])
																	/pow(EvErr_avg_Vn[n], 2), 2);
				err_Vn_6part[bin][n] += pow(Vn_6part[bin][n]
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
			err_Vn_4part[bin][n] = sqrt(err_Vn_4part[bin][n]*(evs - 1)/evs);
      err_Vn_4part_4th[bin][n] = sqrt(err_Vn_4part_4th[bin][n]*(evs - 1)/evs);
			err_Vn_6part[bin][n] = sqrt(err_Vn_6part[bin][n]*(evs - 1)/evs);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    output.push_back({(0.5 + bin)*(100/binned_events.size()),
      Vn_4part[bin][2], err_Vn_4part[bin][2], pow(Vn_4part[bin][2], 4), err_Vn_4part_4th[bin][2], Vn_6part[bin][2], err_Vn_6part[bin][2],
      Vn_4part[bin][3], err_Vn_4part[bin][3], pow(Vn_4part[bin][3], 4), err_Vn_4part_4th[bin][3], Vn_6part[bin][3], err_Vn_6part[bin][3],
      Vn_4part[bin][4], err_Vn_4part[bin][4], pow(Vn_4part[bin][4], 4), err_Vn_4part_4th[bin][4], Vn_6part[bin][4], err_Vn_6part[bin][4],
      Vn_4part[bin][5], err_Vn_4part[bin][5], pow(Vn_4part[bin][5], 4), err_Vn_4part_4th[bin][5], Vn_6part[bin][5], err_Vn_6part[bin][5]
    });
	}

}

vector<vector<double>> Calculator::Calculate_2Particle_Cummulants()
{
  vector<vector<double>> Vn_2part, EvErr_avg_Vn, err_Vn_2part, output;
  Vn_2part.resize(binned_events.size(), vector<double>(7, 0));
  EvErr_avg_Vn.resize(binned_events.size(), vector<double>(7, 0));
  err_Vn_2part.resize(binned_events.size(), vector<double>(7, 0));

	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values

		//	Calculate 2, 4, and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
			Vn_2part[bin][n] = sqrt(avg_Vn[bin][n]);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors

		for(int ev = 0; ev < evs; ev++)
		{
      // Calculate error by event of average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        EvErr_avg_Vn[bin][n] = (avg_Vn[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
      }

			//	Calculate error for 2, 4, and 6 particle cummulants
		 	for (int n = 2; n <= 5; n++)
			{
		 		err_Vn_2part[bin][n] += pow(Vn_2part[bin][n] - sqrt(EvErr_avg_Vn[bin][n]), 2);
			}
		}

		//	Normalize errors
		for (int n = 2; n <= 5; n++)
		{
			err_Vn_2part[bin][n] = sqrt(err_Vn_2part[bin][n]*(evs - 1)/evs);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    output.push_back({(0.5 + bin)*(100/binned_events.size()),
      Vn_2part[bin][2], err_Vn_2part[bin][2],
      Vn_2part[bin][3], err_Vn_2part[bin][3],
      Vn_2part[bin][4], err_Vn_2part[bin][4],
      Vn_2part[bin][5], err_Vn_2part[bin][5]
    });
	}

  return output;
}

vector<vector<double>> Calculator::Calculate_4and6Particle_Cummulants()
{
  vector<vector<double>> Vn_4part, Vn_6part, err_Vn_4part, err_Vn_6part, output;
  Vn_4part.resize(binned_events.size(), vector<double>(7, 0));
  Vn_6part.resize(binned_events.size(), vector<double>(7, 0));
  err_Vn_4part.resize(binned_events.size(), vector<double>(7, 0));
  err_Vn_6part.resize(binned_events.size(), vector<double>(7, 0));

	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values

		//	Calculate 2, 4, and 6 particle cummulants
		for (int n = 2; n <= 5; n++)
		{
			Vn_4part[bin][n] = 2.*pow(avg_Vn[bin][n], 2) - avg_Vn_2nd[bin][n];
			Vn_6part[bin][n] = 0.25*(avg_Vn_3rd[bin][n] - 9.*avg_Vn[bin][n]*avg_Vn_2nd[bin][n] + 12.*pow(avg_Vn[bin][n], 3));
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
    double EvErr_avg_Vn[7] = {0}, EvErr_avg_Vn_2nd[7] = {0},  EvErr_avg_Vn_3rd[7] = {0};

		for(int ev = 0; ev < evs; ev++)
		{
      // Calculate error by event of average Vn, Vn^2, Vn^3
      for (int n = 2; n <= 5; n++)
      {
        EvErr_avg_Vn[n] = (avg_Vn[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 2))/(evs - 1);
        EvErr_avg_Vn_2nd[n] = (avg_Vn_2nd[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 4))/(evs - 1);
        EvErr_avg_Vn_3rd[n] = (avg_Vn_3rd[bin][n]*evs - pow(binned_events[bin][ev].eccentricity[n], 6))/(evs - 1);
      }

			//	Calculate error for 2, 4, and 6 particle cummulants
		 	for (int n = 2; n <= 5; n++)
			{
				err_Vn_4part[bin][n] += pow(Vn_4part[bin][n] - (2.*pow(EvErr_avg_Vn[n], 2) - EvErr_avg_Vn_2nd[n]), 2);
				err_Vn_6part[bin][n] += pow(Vn_6part[bin][n] - 0.25*(EvErr_avg_Vn_3rd[n] - 9.*EvErr_avg_Vn[n]*EvErr_avg_Vn_2nd[n] + 12.*pow(EvErr_avg_Vn[n], 3)), 2);
			}

		}

		//	Normalize errors
		for (int n = 2; n <= 5; n++)
		{
			err_Vn_4part[bin][n] = sqrt(err_Vn_4part[bin][n]*(evs - 1)/evs);
			err_Vn_6part[bin][n] = sqrt(err_Vn_6part[bin][n]*(evs - 1)/evs);
		}
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    output.push_back({(0.5 + bin)*(100/binned_events.size()),
      Vn_4part[bin][2], err_Vn_4part[bin][2], Vn_6part[bin][2], err_Vn_6part[bin][2],
      Vn_4part[bin][3], err_Vn_4part[bin][3], Vn_6part[bin][3], err_Vn_6part[bin][3],
      Vn_4part[bin][4], err_Vn_4part[bin][4], Vn_6part[bin][4], err_Vn_6part[bin][4],
      Vn_4part[bin][5], err_Vn_4part[bin][5], Vn_6part[bin][5], err_Vn_6part[bin][5]
    });
	}

  return output;
}

vector<vector<double>> Calculator::Calculate_V2_V3_CummulantRatio()
{
  vector<double> V2_2part_V3_2part, EvErr_avg_V2, EvErr_avg_V3, err_V2_2part_V3_2part;
  vector<vector<double>> output;
  V2_2part_V3_2part.resize(binned_events.size(), 0);
  EvErr_avg_V2.resize(binned_events.size(), 0);
  EvErr_avg_V3.resize(binned_events.size(), 0);
  err_V2_2part_V3_2part.resize(binned_events.size(), 0);

  for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values

		V2_2part_V3_2part[bin] = sqrt(avg_Vn[bin][2]/avg_Vn[bin][3]);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors

		for(int ev = 0; ev < evs; ev++)
		{
      // Calculate error by event of average Vn, Vn^2, Vn^3
      EvErr_avg_V2[bin] = (avg_Vn[bin][2]*evs - pow(binned_events[bin][ev].eccentricity[2], 2))/(evs - 1);
      EvErr_avg_V3[bin] = (avg_Vn[bin][3]*evs - pow(binned_events[bin][ev].eccentricity[3], 2))/(evs - 1);

			err_V2_2part_V3_2part[bin] += pow(V2_2part_V3_2part[bin] - sqrt(EvErr_avg_V2[bin]/EvErr_avg_V3[bin]), 2);
		}

		//	Normalize errors
		err_V2_2part_V3_2part[bin] = sqrt(err_V2_2part_V3_2part[bin]*(evs - 1)/evs);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    output.push_back({(0.5 + bin)*(100/binned_events.size()), V2_2part_V3_2part[bin], err_V2_2part_V3_2part[bin]});
	}

  return output;
}

vector<vector<double>> Calculator::CalculatePtObservables()
{
  vector<double> mean_pt, mean_pt_2nd, RMS_pt, err_mean_pt, err_mean_pt_2nd, err_RMS_pt;
  vector<vector<double>> output;
  mean_pt.resize(binned_events.size(), 0);
  mean_pt_2nd.resize(binned_events.size(), 0);
  RMS_pt.resize(binned_events.size(), 0);
  err_mean_pt.resize(binned_events.size(), 0);
  err_mean_pt_2nd.resize(binned_events.size(), 0);
  err_RMS_pt.resize(binned_events.size(), 0);

	for (int bin = 0; bin < binned_events.size(); bin++)
	{
		int evs = binned_events[bin].size();

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate values
		for(int ev = 0; ev < evs; ev++) mean_pt[bin] += binned_events[bin][ev].average_pt;
		mean_pt[bin] /= evs;

      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //	Calculate Values
      for(int ev = 0; ev < evs; ev++)
      {
        //	Calculate average Vn, Vn^2, Vn^3
        mean_pt_2nd[bin] += pow(binned_events[bin][ev].average_pt - mean_pt[bin], 2.);
      }


		mean_pt_2nd[bin] /= evs;
		RMS_pt[bin] = sqrt(mean_pt_2nd[bin])/mean_pt[bin];
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//	Calculate errors
		for(int ev = 0; ev < evs; ev++)
		{
			double EvErr_mean_pt = (mean_pt[bin]*evs - binned_events[bin][ev].average_pt)/(evs - 1.);
			double EvErr_mean_pt_2nd = (mean_pt_2nd[bin]*evs - pow(binned_events[bin][ev].average_pt - EvErr_mean_pt, 2))/(evs - 1);

		 	err_mean_pt[bin] += pow(mean_pt[bin] - EvErr_mean_pt, 2);
		 	err_mean_pt_2nd[bin] += pow(mean_pt_2nd[bin] - EvErr_mean_pt_2nd, 2);
		 	err_RMS_pt[bin] = pow(RMS_pt[bin] - sqrt(EvErr_mean_pt_2nd)/EvErr_mean_pt, 2);
		}

		err_mean_pt[bin] = sqrt(err_mean_pt[bin]*(evs - 1)/evs);
		err_mean_pt_2nd[bin] = sqrt(err_mean_pt_2nd[bin]*(evs - 1)/evs);
		err_RMS_pt[bin] = sqrt(err_RMS_pt[bin]*(evs - 1)/evs);
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    output.push_back({(0.5 + bin)*(100/binned_events.size()),
      mean_pt[bin], err_mean_pt[bin],
      mean_pt_2nd[bin], err_mean_pt_2nd[bin],
      RMS_pt[bin], err_RMS_pt[bin]
    });
	}

  return output;
}
