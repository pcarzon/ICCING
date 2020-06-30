#include "event.h"

Event::Event()
{

}

Event::~Event()
{

}

Event::Event(const Event &original)
{

}

void Event::CopyEvent(const Event &e)
{

}

Event& Event::operator= (const Event& original)
{
	CopyEvent(original);
	return *this;
}

void Event::SampleEnergy(int location, int radius)
{

}

void Event::UpdateDensity(int densityType, vector<vector<double>> inputDensity)
{

}
