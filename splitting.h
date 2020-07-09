#ifndef Splitter_H
#define Splitter_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "global.h"

using namespace std;

class Splitter
{
private:

	void CopySplitter(const Splitter &e);

public:

  Splitter();
  ~Splitter();

  Splitter(const Splitter &original);
  Splitter& operator=(const Splitter& original);

  Quarks SplitSample(Sample sampled_energy);  //  See: SplittingSample in ICCING_v0_1_8.nb
};
#endif
