#include "ecc.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Eccentricity
//##########################################################################################
Eccentricity::Eccentricity()
{
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
Eccentricity::~Eccentricity()
{
  CleanEccentricity();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Eccentricity::Eccentricity(const Eccentricity &original)
{
  CopyEccentricity(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Eccentricity Copy Function
//##########################################################################################
void Eccentricity::CopyEccentricity(const Eccentricity &e)
{
  x_center_of_mass = e.x_center_of_mass;
  y_center_of_mass = e.y_center_of_mass;
  sparse_density = e.sparse_density;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Eccentricity& Eccentricity::operator= (const Eccentricity& original)
{
	CopyEccentricity(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Calculate eccentricity
//##########################################################################################
vector<double> Eccentricity::StandardCalculation(string density_type, int m, int n)
{
  int column;
  int max = sparse_density.size();
  vector <double> distance_squared, phi;
  double psi, radius;
	double psi_top = 0, psi_bottom = 0, normalization = 0, x_component = 0, y_component = 0, weight = 0;
  double eccentricity = 0, etot = 0;

  distance_squared.resize(max);
  phi.resize(max);

  if (density_type == "Energy") { column = 2; }

	for (int s=0;s<max;s++)
  {
	   x_component = (sparse_density[s][0] - x_center_of_mass);
	   y_component = (sparse_density[s][1] - y_center_of_mass);
	   distance_squared[s] = pow(x_component, 2) + pow(y_component, 2);

     weight = sparse_density[s][column]*pow(distance_squared[s], (m/2.));
     normalization += weight;

     phi[s] = atan2(y_component, x_component); // angle of fluid cells

     psi_top += weight*sin(1.0*n*phi[s]);
	   psi_bottom += weight*cos(1.0*n*phi[s]);

     etot += sparse_density[s][column];
	}
  // m is radial weight
  // n is anglular weight

	psi_top /= max;
	psi_bottom /= max;

  // relative event plane angle (perp to major axis) (coming out of flat sides of shape)
	psi = 1./(1.0*n)*atan2(psi_top, psi_bottom);

	for (int s=0;s<max;s++)
  {
    eccentricity += sparse_density[s][column]*pow(distance_squared[s], m/2.)*cos(n*(phi[s] - psi));
  }

  eccentricity /= normalization;

  // top and bottom of eccentricity is technically divided by number of particles (max)
	radius = normalization/etot;

	return {eccentricity, psi, radius};

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Calculate eccentricities, seperating positive and negative density values
//##########################################################################################
vector<double> Eccentricity::NewCalculation(string density_type, int m, int n)
{
  int column;
  int max = sparse_density.size(), max_pos = 0, max_neg = 0;
  vector <double> distance_squared, phi;
  double psi_pos, psi_neg, radius_pos, radius_neg;
	double psi_top_pos = 0, psi_bottom_pos = 0, psi_top_neg = 0, psi_bottom_neg = 0, normalization_pos = 0, normalization_neg = 0;
  double x_component = 0, y_component = 0, weight = 0;
  double eccentricity_pos = 0, etot_pos = 0, eccentricity_neg = 0, etot_neg = 0;

  distance_squared.resize(max);
  phi.resize(max);

  if (density_type == "Baryon") { column = 3; }
  else if (density_type == "Strange") { column = 4; }
  else if (density_type == "Charge") { column = 5; }

	for (int s=0;s<max;s++)
  {
	   x_component = (sparse_density[s][0] - x_center_of_mass);
	   y_component = (sparse_density[s][1] - y_center_of_mass);
	   distance_squared[s] = pow(x_component, 2) + pow(y_component, 2);

     weight = sparse_density[s][column]*pow(distance_squared[s], (m/2.));

     if (sparse_density[s][column] < 0)
     {  normalization_neg += weight;  }
     else if (sparse_density[s][column] > 0)
     {  normalization_pos += weight;  }

     phi[s] = atan2(y_component, x_component); // angle of fluid cells

     if (sparse_density[s][column] < 0)
     {
       psi_top_neg += weight*sin(1.0*n*phi[s]);
     	 psi_bottom_neg += weight*cos(1.0*n*phi[s]);
       max_neg++;
     }
     else if (sparse_density[s][column] > 0)
     {
       psi_top_pos += weight*sin(1.0*n*phi[s]);
     	 psi_bottom_pos += weight*cos(1.0*n*phi[s]);
       max_pos++;
     }

     if (sparse_density[s][column] < 0)
     {  etot_neg += sparse_density[s][column];  }
     else if (sparse_density[s][column] > 0)
     {  etot_pos += sparse_density[s][column];  }
	}
  // m is radial weight
  // n is anglular weight

    psi_top_neg /= max_neg;
    psi_bottom_neg /= max_neg;

    psi_top_pos /= max_pos;
    psi_bottom_pos /= max_pos;

  // relative event plane angle (perp to major axis) (coming out of flat sides of shape)
  // need to divide parts by normalization so atan2 gets the quadrant correct and thus give positive eccs
	psi_neg = 1./(1.0*n)*atan2(psi_top_neg/normalization_neg, psi_bottom_neg/normalization_neg);

  psi_pos = 1./(1.0*n)*atan2(psi_top_pos/normalization_pos, psi_bottom_pos/normalization_pos);

	for (int s=0;s<max;s++)
  {
    if (sparse_density[s][column] < 0)
    {  eccentricity_neg += sparse_density[s][column]*pow(distance_squared[s], m/2.)*cos(n*(phi[s] - psi_neg));  }
    else if (sparse_density[s][column] > 0)
    {  eccentricity_pos += sparse_density[s][column]*pow(distance_squared[s], m/2.)*cos(n*(phi[s] - psi_pos));  }
  }

  eccentricity_neg /= normalization_neg;
  eccentricity_pos /= normalization_pos;

  // top and bottom of eccentricity is technically divided by number of particles (max)
	radius_neg = normalization_neg/etot_neg;
  radius_pos = normalization_pos/etot_pos;

	return {eccentricity_neg, psi_neg, radius_neg, eccentricity_pos, psi_pos, radius_pos};
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Calculate All eccentricities for given event
//##########################################################################################
vector<vector<vector<double>>> Eccentricity::CalculateEccentricities(int grid_max, double grid_step, vector<vector<vector<double>>> density)
{
  double x, y, energy = 0;

  //******************************************************************************************
  //  Take full density grid and convert to sparse density structure for easy and quick processing
  //******************************************************************************************
  for (int i = 0; i < density[0].size(); i++)
  {
    for (int j = 0; j < density[0][0].size(); j++)
    {
      if (density[0][i][j] != 0)
      {
        x = -grid_max + i*grid_step;  //  Converts grid point to physical x-value
        y = -grid_max + j*grid_step;  //  Converts grid point to physical y-value

        x_center_of_mass += x*density[0][i][j];
        y_center_of_mass += y*density[0][i][j];
        energy += density[0][i][j];
        sparse_density.push_back({x, y, density[0][i][j], density[1][i][j], density[2][i][j], density[3][i][j]});
      }
    }
  }

  x_center_of_mass /= energy;
  y_center_of_mass /= energy;

  //******************************************************************************************
  //  Calculate eccentricities and return in structure for easy output
  //******************************************************************************************
  return {{StandardCalculation("Energy",2,2), StandardCalculation("Energy",3,3), StandardCalculation("Energy",4,4), StandardCalculation("Energy",5,5)}
         ,{NewCalculation("Baryon",2,2), NewCalculation("Baryon",3,3), NewCalculation("Baryon",4,4), NewCalculation("Baryon",5,5)}
         ,{NewCalculation("Strange",2,2), NewCalculation("Strange",3,3), NewCalculation("Strange",4,4), NewCalculation("Strange",5,5)}
         ,{NewCalculation("Charge",2,2), NewCalculation("Charge",3,3), NewCalculation("Charge",4,4), NewCalculation("Charge",5,5)}};
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Clean class
//##########################################################################################
void Eccentricity::CleanEccentricity()
{
  x_center_of_mass = 0;
  y_center_of_mass = 0;
  sparse_density.clear();
}
//__________________________________________________________________________________________
