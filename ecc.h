#ifndef Eccentricity_H
#define Eccentricity_H

#include <Eccentricitystream>
#include <string>
#include <cmath>
#include <vector>

#include "event.h"

using namespace std;

class Eccentricity
{
private:

	void CopyEccentricity(const Eccentricity &e);

public:

  Eccentricity();
  ~Eccentricity();

  Eccentricity(const Eccentricity &original);
  Eccentricity& operator=(const Eccentricity& original);

};
#endif
