#include <iostream>

using namespace std;


/*
�����������ķ���
*/

double Norm(double *a, double *b, int n);
/*
Jacobi��������ʵ����AX=b
*/

double* Jacobi(int n, double *x0,double **a, int N, double TOL);
void testJacobi();

/*
Gauss-Seidel��������ʵ����AX=b
*/

double* Gauss_Seidel(int n, double *x0,double **a, int N, double TOL);
void testGauss_Seidel();

/*
SOR�㷨
*/

double* SOR(double w,int n, double *x0,double **a, int N, double TOL);
void testSOR();