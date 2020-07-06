#ifndef Eccentricity_H
#define Eccentricity_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

class Event;

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
