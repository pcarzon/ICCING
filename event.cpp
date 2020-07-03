#include "event.h"

Event::Event()
{
  initial_energy.resize(401, vector<double>(401, 0));
}

Event::~Event()
{
  initial_energy.clear();
  density.clear();
  sampled_energy.clear();

  total_initial_energy = 0;
  total_energy = 0;
  seed = 0;
}

Event::Event(const Event &original)
{
	CopyEvent(original);
}

void Event::CopyEvent(const Event &e)
{
  initial_energy = e.initial_energy;
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
  initial_energy = initEnergy;
}

void Event::SampleEnergy(int location, double radius)
{

}

void Event::UpdateDensity(int densityType, vector<vector<double>> inputDensity)
{

}
