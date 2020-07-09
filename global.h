#ifndef Global_H
#define Global_H

#include<vector>

using namespace std;

struct Charge
{
private:
  vector<double> gluon = {0., 0., 0., 0.};
  vector<double> up = {0.0023, 0.33333333, 0, 0.66666666};
  vector<double> down = {0.0048, 0.33333333, 0, -0.33333333};
  vector<double> strange = {0.095, 0.33333333, -1, -0.33333333};
  vector<double> charm = {1.29, 0.33333333, 0, 0.66666666};

  vector<double> current_charge;
public:
  void Gluon() {  current_charge = gluon; }
  void Up() { current_charge = up;  }
  void Down() { current_charge = down;  }
  void Strange() { current_charge = strange;  }
  void Charm() { current_charge = charm;  }

  vector<double> GetCharge() {  return current_charge;  }
}

struct Quarks
{
private:
  Charge flavor;

  double alpha;
  vector<double> position[2];

public:
  void CreateQuarks (Charge flavor_, double alpha_, double delta_x_, double delta_y_)
  {
    flavor = flavor_; alpha = alpha_; position[0] = delta_x_; position[1] = delta_y_;
  }
  double GetAlpha () {  return alpha; }
  vector<double> GetPosition () {  return position; }
}

struct Sample
{
  double e_tot;
  double q_s;
}

#endif
