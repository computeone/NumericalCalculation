#include <iostream>
#include <vector>

using namespace std;

/*
积分的数值算法
*/
/*
复合Simpson法则
*/
double CompoundSimpson(double a, double b, int n);
void testCompoundSimpson();
double function(double d);

/*
Romberg外推法实现
*/

double** Romberg(double a, double b, int n,double(*function)(double));
void testRomberg();
double functionRomberg(double x);

/*
自适应求积分
*/

double AdaptiveIntegration(double a, double b, double tol, int N, double(*function)(double));
void testAdaptiveIntegration();
double functionAdaptive(double x);