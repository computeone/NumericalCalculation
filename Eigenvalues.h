#include <iostream>

using namespace std;

/*
�ƽ�����ֵ
*/

/*
�ݷ�
*/

double PowerMethod(double *x0, double **a,int n, double TOL, int N);
void testPowerMethod();

/*
�Գ��ݷ�
*/

double SymmetryPowerMethod(double *x0, double **a, int n, double TOL, int N);
void testSymmetryPowerMethod();

/*
���ݷ�
*/
double AntiPowerMethod(double *x0, double **a, int n, double TOL, int N);
void testAntiPowerMethod();