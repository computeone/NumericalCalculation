#include <iostream>

using namespace std;;


//Horner算法，求多项式p(x)=AnX^n+An-1X^n-1+....a1X+a0=(x-x0)Q(x)+b0的值和导数

pair<double,double> Horner(double x0);
pair<double,double> Horner(double x0, int n,double a[],double(*function)(double));

double testMuller(double x);
// Muller算法，求多项式的m重根
double Muller(double x0, double x1, double x2, double tol, int n, double(*function)(double));
