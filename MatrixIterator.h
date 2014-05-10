#include <iostream>

using namespace std;


/*
求两个向量的范数
*/

double Norm(double *a, double *b, int n);
/*
Jacobi迭代法，实现求AX=b
*/

double* Jacobi(int n, double *x0,double **a, int N, double TOL);
void testJacobi();

/*
Gauss-Seidel迭代法，实现求AX=b
*/

double* Gauss_Seidel(int n, double *x0,double **a, int N, double TOL);
void testGauss_Seidel();

/*
SOR算法
*/

double* SOR(double w,int n, double *x0,double **a, int N, double TOL);
void testSOR();