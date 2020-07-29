#ifndef Functions_H
#define Functions_H

#include <iostream>
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

//    friend bool operator< (const SplineSet& lhs, const double& rhs)
//    {
//      return lhs.x < rhs;
//    }
};
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Natural Cubic Spline
//	c++ implimentation: https://bit.ly/3hP8dUd
//	algorithm: https://bit.ly/3gdXZfK
//##########################################################################################
vector<SplineSet> CubicSpline(vector<double> &x, vector<double> &y);

double InterpolateValue(SplineSet range, double value);

SplineSet FindRange(vector<SplineSet> function, double value);

#endif
