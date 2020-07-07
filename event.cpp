#include "event.h"

Event::Event()
{

}

Event::~Event()
{
  CleanEvent();
}

Event::Event(const Event &original)
{
	CopyEvent(original);
}

void Event::CopyEvent(const Event &e)
{
  initial_energy = e.initial_energy;
  t_a = e.t_a;
  t_b = e.t_b;
  density = e.density;
  sampled_energy = e.sampled_energy;

  total_initial_energy = e.total_initial_energy;
  total_energy = e.total_energy;
  seed = e.seed;
}

Event& Event::operator= (const Event& original)
{
	CopyEvent(original);
	return *this;
}

void Event::ReadInitialEnergy(vector<vector<double>> initEnergy)
{
  initial_energy.resize(401, vector<double>(401, 0));
  initial_energy = initEnergy;
}

double Event::SampleEnergy(int location, double radius)
{

}

void Event::UpdateDensity(int densityType, vector<vector<double>> inputDensity)
{

}

void Event::CleanEvent()
{
  initial_energy.clear();
  t_a.clear();
  t_b.clear();
  density.clear();
  sampled_energy.clear();
  /*vector<vector<double>>().swap(initial_energy);
  vector<vector<double>>().swap(t_a);
  vector<vector<double>>().swap(t_b);
  vector<vector<vector<double>>>().swap(density);
  vector<vector<double>>().swap(sampled_energy);*/

  total_initial_energy = 0;
  total_energy = 0;
  seed = 0;
}
