#include "ecc.h"

Eccentricity::Eccentricity()
{

}

Eccentricity::~Eccentricity()
{

}

Eccentricity::Eccentricity(const Eccentricity &original)
{
  CopyEccentricity(original);
}

void Eccentricity::CopyEccentricity(const Eccentricity &e)
{

}

Eccentricity& Eccentricity::operator= (const Eccentricity& original)
{
	CopyEccentricity(original);
	return *this;
}
