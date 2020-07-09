#ifndef ProbDist_H
#define ProbDist_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

class ProbDist
{
private:

	void CopyProbDist(const ProbDist &e);

public:

  ProbDist();
  ~ProbDist();

  ProbDist(const ProbDist &original);
  ProbDist& operator=(const ProbDist& original);

};
#endif
