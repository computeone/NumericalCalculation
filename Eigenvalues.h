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

/*
Wielandt收缩法
*/

void Wielandt(double *v, double **a, int n, double r0, double TOL, int N);
void testWielandt();

/*
Householder变换算法
*/
void Householder(double **a, int n);
void testHouseholder();


/*
QR分解算法
*/

double QR(double **a, double *b, int n,double TOL,int M);
void testQR();
