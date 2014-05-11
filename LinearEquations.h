#include <iostream>
#include <iomanip>
using namespace std;

/*
向后代换的Guass消元法
*/

double* GaussianElimination(double **A, int n);
void testGaussianElimination();

/*
部分选主元Gauss消去法
*/

double* GaussianPivotElimination(double **a, int n);
void testGaussianPivotElimination();
/*
LU分解算法
*/

void LU(int n, double **a);
void testLU();