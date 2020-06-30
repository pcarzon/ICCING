#ifndef IO_H
#define IO_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "event.h"

using namespace std;

class IO
{
private:

	void CopyIO(const IO &e);

public:
  
  IO();
  ~IO();

  IO(const IO &original);
  IO& operator=(const IO& original);

};
#endif
