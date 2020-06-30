#ifndef IO_H
#define IO_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "event.h"
#include "ecc.h"

using namespace std;

class IO
{
private:

  string input_file;
  string output_file;
  int input_type;
  int output_type;

	void CopyIO(const IO &e);

  void OutputFullDensityGrids(const Event &event);
  void OutputSparseDensityGrids(const Event &event);
  void OutputEccentricities(const Eccentricity &ecc);

public:

  IO(string input, string output, int iType, int oType);
  ~IO();

  IO(const IO &original);
  IO& operator=(const IO& original);

  Event ReadEvent();
  void WriteEvent(const Event &event);
};
#endif
