#include <iostream>
#include <iomanip>
using namespace std;

/*
��������Guass��Ԫ��
*/

double* GaussianElimination(double **A, int n, int k);
void testGaussianElimination();

/*
LU�ֽ��㷨
*/

void LU(int n, double **a);
void testLU();