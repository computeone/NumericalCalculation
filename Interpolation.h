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
/*
Hermite多项式插值
*/

double** Hermite(vector<double> xn, vector<double> fx,vector<double> fxx);

/*
自然三次样条插值
*/

double** CubicSpline(vector<double> xn, vector<double> fx);

//用数据测试三次样条插值
void testCubicSpline();