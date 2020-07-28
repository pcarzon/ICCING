#include "functions.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Natural Cubic Spline
//	c++ implimentation: https://bit.ly/3hP8dUd
//	algorithm: https://bit.ly/3gdXZfK
//##########################################################################################
vector<SplineSet> CubicSpline(vector<double> &x, vector<double> &y)
{
    int n = x.size()-1;
    vector<double> a;
    a.insert(a.begin(), y.begin(), y.end());
    vector<double> b(n);
    vector<double> d(n);
    vector<double> h;

    for(int i = 0; i < n; ++i)
        h.push_back(x[i+1]-x[i]);

    vector<double> alpha;
    alpha.push_back(0);
    for(int i = 1; i < n; ++i)
        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

    vector<double> c(n+1);
    vector<double> l(n+1);
    vector<double> mu(n+1);
    vector<double> z(n+1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for(int i = 1; i < n; ++i)
    {
        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for(int j = n-1; j >= 0; --j)
    {
        c[j] = z [j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
        d[j] = (c[j+1]-c[j])/3/h[j];
    }

    vector<SplineSet> output_set(n);
    for(int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }
    return output_set;
}

SplineSet SplineInterval(double value)
{

}

double InterpolateValue(SplineSet poly, double value)
{
  double change_x = value - poly.x;
  return (poly.a + poly.b*change_x + poly.c*pow(change_x, 2) + poly.d*pow(change_x, 3));
}
