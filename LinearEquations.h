#include <iostream>
#include <iomanip>
using namespace std;

/*
��������Guass��Ԫ��
*/

double* GaussianElimination(double **A, int n);
void testGaussianElimination();

/*
����ѡ��ԪGauss��ȥ��
*/

double* GaussianPivotElimination(double **a, int n);
void testGaussianPivotElimination();
/*
LU�ֽ��㷨
*/

void LU(int n, double **a);
void testLU();