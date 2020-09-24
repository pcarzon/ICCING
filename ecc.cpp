#include "ecc.h"

vector<double> Eccentricities(vector<vector<double>> grid, double grid_step)
{
  cout << "started eccs" << endl;

  double x_center_of_mass = 0, y_center_of_mass = 0, total = 0;
  complex<double> epsilon_2 = {0, 0}, epsilon_3 = {0, 0}, epsilon_4 = {0, 0}, epsilon_5 = {0, 0};

  complex<double> numerator2 = {0, 0}, denominator2 = {0, 0};
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

      numerator2 += pow(complex<double>(del_x, del_y), 2.)*grid[p][q];
      denominator2 += pow(pow(del_x, 2.) + pow(del_y, 2.), 2/2)*grid[p][q];

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
  cout << numerator2 << " " << denominator2 << " " << numerator3 << " " << denominator3 << " " << numerator4 << " " << denominator4 << endl;
  epsilon_2 = -numerator2/denominator2;
  epsilon_3 = -numerator3/denominator3;
  epsilon_4 = -numerator4/denominator4;
  epsilon_5 = -numerator5/denominator5;
  cout << "got eccs" << endl;
  cout << abs(epsilon_2) << " " << abs(epsilon_3) << " " << abs(epsilon_4) << " " << abs(epsilon_5) << endl;

  return {abs(epsilon_2), abs(epsilon_3), abs(epsilon_4), abs(epsilon_5)};
}

vector<double> NewEccentricities(vector<vector<double>> grid, double grid_step)
{

}
