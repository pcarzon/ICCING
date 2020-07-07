#ifndef Splitter_H
#define Splitter_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "event.h"

using namespace std;

class Event;

class Splitter
{
private:

	void CopySplitter(const Splitter &e);

public:

  Splitter();
  ~Splitter();

  Splitter(const Splitter &original);
  Splitter& operator=(const Splitter& original);

  void Split(double q_s);
};
#endif
