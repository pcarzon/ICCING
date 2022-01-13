#ifndef Global_H
#define Global_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <vector>
#include <random>

using namespace std;


//__________________________________________________________________________________________
//##########################################################################################
//  Data Structure Used to keep track of particles and their charges
//##########################################################################################
struct Charge
{
private:
  //  Predefined mass and charge vectors for gluon and quarks
  //    particle_type[BSQ=0,UDS=1][Mass=0,baryon#=1,strangeness=2,em_charge=3]
  vector<vector<double>> gluon = {{0., 0., 0., 0.},{0., 0., 0., 0.}};
  vector<vector<double>> up = {{0.0023, 0.3333, 0., 0.6666},{0.0023, 1., 0., 0.}};
  vector<vector<double>> down = {{0.0048, 0.3333, 0., -0.3333},{0.0048, 0., 1., 0.}};
  vector<vector<double>> strange = {{0.095, 0.3333, -1., -0.3333},{0.095, 0., 0., 1.}};
  vector<vector<double>> charm = {{1.29, 0.3333, 0., 0.6666},{1.29, 0., 0., 0.}};

  //  Stores the currently set particle
  vector<double> current_charge;

public:
  //  Function to set the charge of the Charge object
  void Gluon(int charge_type) {  current_charge = gluon[charge_type]; }
  void Up(int charge_type) { current_charge = up[charge_type];  }
  void Down(int charge_type) { current_charge = down[charge_type];  }
  void Strange(int charge_type) { current_charge = strange[charge_type];  }
  void Charm(int charge_type) { current_charge = charm[charge_type];  }

  //  Returns the set charge of the object
  vector<double> GetCharge() {  return current_charge;  }
};
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Data Structure for generated quark pairs
//##########################################################################################
struct Quarks
{
private:
  Charge flavor;  //  Stores what the quark charges are

  double energy_fraction; //  Energy fraction of quarks
  double alpha; //  Momentum fraction of quark (fraction for anit-quark = 1 - alpha)
  vector<double> position{0, 0};  //  Stores the position of the quark pair

public:
  //  Create quarks with given flavor, momentum fraction, and position
  void CreateQuarks (Charge flavor_, double energy_fraction_, double alpha_, double delta_x_, double delta_y_)
  {
    flavor = flavor_; energy_fraction = energy_fraction_; alpha = alpha_; position[0] = delta_x_; position[1] = delta_y_;
  }
  double GetEnergyFraction () {  return energy_fraction; } //  Returns energy fraction
  double GetAlpha () {  return alpha; } //  Returns momentum fraction
  vector<double> GetPosition () {  return position; } //  Returns quark pair position
  vector<double> GetCharge() {  return flavor.GetCharge();  }
};
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Data Structure for sample of initial energy density
//##########################################################################################
struct Sample
{
  double e_tot;
  double q_s;
};
//__________________________________________________________________________________________

#endif
