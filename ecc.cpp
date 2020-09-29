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
//     if (isnan(normalization))
//     cout << "y1_component " << sparse_density[1][s] << " y2_component " << y_center_of_mass << endl;
//  if (sparse_density[s][column] > 24) cout << "y " << y << " sparse y " << sparse_density[sparse_density.size()][1] << endl;

     phi[s] = atan2(y_component, x_component); // angle of fluid cells

     psi_top += weight*sin(1.0*n*phi[s]);
	   psi_bottom += weight*cos(1.0*n*phi[s]);

     etot += sparse_density[s][column];
	}
//  cout << phi[1] << " " << phi[2] << " " << phi[3] << endl;
//  cout << normalization << " " << psi_top << " " << psi_bottom << endl;
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
//  cout << "eccentricity " << eccentricity << endl;
  eccentricity /= normalization;
//  cout << "eccentricity/normalized " << eccentricity << endl;

  // top and bottom of eccentricity is technically divided by number of particles (max)
	radius = normalization/etot;

	return {eccentricity, psi, radius};

}
/*
vector<double> Eccentricities(vector<vector<double>> grid, double grid_step)
{
  cout << "started eccs" << endl;

  double x_center_of_mass = 0, y_center_of_mass = 0, total = 0;
  complex<double> epsilon_2 = {0, 0}, epsilon_3 = {0, 0}, epsilon_4 = {0, 0}, epsilon_5 = {0, 0};

  complex<double> numeratodistance_squared = {0, 0}, denominatodistance_squared = {0, 0};
  complex<double> numerator3 = {0, 0}, denominator3 = {0, 0};
  complex<double> numerator4 = {0, 0}, denominator4 = {0, 0};
  complex<double> numerator5 = {0, 0}, denominator5 = {0, 0};

  for (int p = 0; p < grid.size(); p++)
  {
    for (int q = 0; q < grid[0].size(); q++)
    {
      x_center_of_mass += p*grid[p][q];
      y_center_of_mass += q*grid[p][q];
      total += grid[p][q];
    }
  }

  x_center_of_mass = (1/total)*x_center_of_mass;
  y_center_of_mass = (1/total)*y_center_of_mass;
  cout << "got center of mass x=" << x_center_of_mass << " y=" << y_center_of_mass << endl;

  double del_x = 0, del_y = 0;
  for (int p = 0; p < grid.size(); p++)
  {
    for (int q = 0; q < grid[0].size(); q++)
    {
      del_x = (p - x_center_of_mass)*grid_step;
      del_y = (q - y_center_of_mass)*grid_step;

      numeratodistance_squared += pow(complex<double>(del_x, del_y), 2.)*grid[p][q];
      denominatodistance_squared += pow(pow(del_x, 2.) + pow(del_y, 2.), 2/2)*grid[p][q];

      numerator3 += pow(complex<double>(del_x, del_y), 3.)*grid[p][q];
//      cout << pow(complex<double>(del_x, del_y), 3.) << endl;
      denominator3 += pow(pow(del_x, 2.) + pow(del_y, 2.), 3/2)*grid[p][q];
//      cout << pow(pow(del_x, 2.) + pow(del_y, 2.), 3/2) << endl;

      numerator4 += pow(complex<double>(del_x, del_y), 4.)*grid[p][q];
      denominator4 += pow(pow(del_x, 2.) + pow(del_y, 2.), 4/2)*grid[p][q];

      numerator5 += pow(complex<double>(del_x, del_y), 5.)*grid[p][q];
      denominator5 += pow(pow(del_x, 2.) + pow(del_y, 2.), 5/2)*grid[p][q];
    }
  }
  cout << "got num and denom" << endl;
  cout << numeratodistance_squared << " " << denominatodistance_squared << " " << numerator3 << " " << denominator3 << " " << numerator4 << " " << denominator4 << endl;
  epsilon_2 = -numeratodistance_squared/denominatodistance_squared;
  epsilon_3 = -numerator3/denominator3;
  epsilon_4 = -numerator4/denominator4;
  epsilon_5 = -numerator5/denominator5;
  cout << "got eccs" << endl;
  cout << abs(epsilon_2) << " " << abs(epsilon_3) << " " << abs(epsilon_4) << " " << abs(epsilon_5) << endl;

  return {abs(epsilon_2), abs(epsilon_3), abs(epsilon_4), abs(epsilon_5)};
}
*/
vector<double> NewEccentricities(vector<vector<double>> grid, double grid_step)
{

}

void Eccentricity::CalculateEccentricities(int grid_max, double grid_step, vector<vector<vector<double>>> density)
{
  double x, y, energy = 0;
  for (int i = 0; i < density[0].size(); i++)
  {
    for (int j = 0; j < density[0][0].size(); j++)
    {
      if (density[0][i][j] != 0)
      {
        x = -grid_max + i*grid_step;  //  Converts grid point to physical x-value
        y = -grid_max + j*grid_step;  //  Converts grid point to physical y-value
      //  if (y > 1) cout << "grid_max " << grid_max << " j " << j << " grid_step " << grid_step << endl;
        x_center_of_mass += x*density[0][i][j];
        y_center_of_mass += y*density[0][i][j];
        energy += density[0][i][j];
        sparse_density.push_back({x, y, density[0][i][j], density[1][i][j], density[2][i][j], density[3][i][j]});
//        if (sparse_density[sparse_density.size()-1][1] > 24) cout << "y " << y << " sparse y " << sparse_density[sparse_density.size()][1] << endl;

      }
    }
  }

  x_center_of_mass /= energy;
  y_center_of_mass /= energy;
  cout << x_center_of_mass << " " << y_center_of_mass << endl;
  cout << StandardCalculation("Energy",2,2)[0] << " " << StandardCalculation("Energy",3,3)[0] << endl;

//  return {StandardCalculation("Energy",2,2), StandardCalculation("Energy",3,3), StandardCalculation("Energy",4,4), StandardCalculation("Energy",5,5)};


//  eccentricities = Eccentricities(density[0], grid_step);

}
