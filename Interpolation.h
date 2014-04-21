#include <iostream>
#include <vector>
using namespace std;

/*
插值与多项式逼近
*/

double** Neville(double x,vector<double> xn, vector<double> fx);

/*
牛顿差商插值法
*/

double** divideNewton(vector<double> xn, vector<double> fx);