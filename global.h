#ifndef Global_H
#define Global_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include<vector>

using namespace std;

//__________________________________________________________________________________________
//##########################################################################################
//  Data Structure Used to keep track of particles and their charges
//##########################################################################################
struct Charge
{
private:
  // BSQ/UDS flag to change charges
  //  Predefined mass and charge vectors for gluon and quarks
  vector<double> gluon = {0., 0., 0., 0.};
  vector<double> up = {0.0023, 0.33333333, 0., 0.66666666};
  vector<double> down = {0.0048, 0.33333333, 0., -0.33333333};
  vector<double> strange = {0.095, 0.33333333, -1, -0.33333333};
  vector<double> charm = {1.29, 0.33333333, 0., 0.66666666};

  //  Stores the currently set particle
  vector<double> current_charge;
public:
  //  Function to set the charge of the Charge object
  void Gluon() {  current_charge = gluon; }
  void Up() { current_charge = up;  }
  void Down() { current_charge = down;  }
  void Strange() { current_charge = strange;  }
  void Charm() { current_charge = charm;  }

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

  double alpha; //  Momentum fraction of quark (fraction for anit-quark = 1 - alpha)
  vector<double> position{0, 0};  //  Stores the position of the quark pair

public:
  //  Create quarks with given flavor, momentum fraction, and position
  void CreateQuarks (Charge flavor_, double alpha_, double delta_x_, double delta_y_)
  {
    flavor = flavor_; alpha = alpha_; position[0] = delta_x_; position[1] = delta_y_;
  }
  double GetAlpha () {  return alpha; } //  Returns momentum fraction
  vector<double> GetPosition () {  return position; } //  Returns quark pair position
  vector<double> GetCharge() {  return flavor.GetCharge;  }
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
