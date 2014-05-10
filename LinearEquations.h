#include <iostream>
#include <iomanip>
using namespace std;

/*
向后代换的Guass消元法
*/

double* GaussianElimination(double **A, int n, int k);
void testGaussianElimination();

/*
LU分解算法
*/

void LU(int n, double **a);
void testLU();