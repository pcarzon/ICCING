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

}

Splitter& Splitter::operator= (const Splitter& original)
{
	CopySplitter(original);
	return *this;
}

Charges Splitter::Split(double q_s)
{

}
