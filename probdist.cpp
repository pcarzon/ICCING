#include "probdist.h"

ProbDist::ProbDist()
{

}

ProbDist::~ProbDist()
{

}

ProbDist::ProbDist(const ProbDist &original)
{
  CopyProbDist(original);
}

void ProbDist::CopyProbDist(const ProbDist &e)
{

}

ProbDist& ProbDist::operator= (const ProbDist& original)
{
	CopyProbDist(original);
	return *this;
}
