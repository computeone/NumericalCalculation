#include <iostream>

using namespace std;

/*
逼近特征值
*/

/*
幂法
*/

double PowerMethod(double *x0, double **a,int n, double TOL, int N);
void testPowerMethod();

/*
对称幂法
*/

double SymmetryPowerMethod(double *x0, double **a, int n, double TOL, int N);
void testSymmetryPowerMethod();

/*
反幂法
*/
double AntiPowerMethod(double *x0, double **a, int n, double TOL, int N);
void testAntiPowerMethod();