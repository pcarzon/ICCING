#ifndef Functions_H
#define Functions_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

//__________________________________________________________________________________________
//##########################################################################################
//  Data Structure for spline functions on interval x - x+1
//##########################################################################################
struct SplineSet{
    double a;
    double b;
    double c;
    double d;
    double x;

    SplineSet& operator=(const SplineSet& s)
    {
      a = s.a;
      b = s.b;
      c = s.c;
      d = s.d;
      x = s.x;
    }

};
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Natural Cubic Spline
//	c++ implimentation: https://bit.ly/3hP8dUd
//	algorithm: https://bit.ly/3gdXZfK
//##########################################################################################
vector<SplineSet> CubicSpline(vector<double> &x, vector<double> &y);

//  Interpolate Value in given SplineSet range
double InterpolateValue(SplineSet range, double value);

//  Find the range that contains the value
SplineSet FindRange(vector<SplineSet> function, double value);

#endif
