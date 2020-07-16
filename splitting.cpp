#include "splitting.h"

Splitter::Splitter()
{

}

Splitter::~Splitter()
{

}

Splitter::Splitter(const Splitter &original)
{
  CopySplitter(original);
}

void Splitter::CopySplitter(const Splitter &e)
{
  flavor_chemistry = e.flavor_chemistry;
	dipole_model = e.dipole_model;
	alpha_s = e.alpha_s;
	alpha_min = e.alpha_min;
	r_max = e.r_max;

}

Splitter& Splitter::operator= (const Splitter& original)
{
	CopySplitter(original);
	return *this;
}

Quarks Splitter::SplitSample(Sample sampled_energy)
{

}
