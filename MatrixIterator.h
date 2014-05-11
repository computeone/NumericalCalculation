#include <iostream>

using namespace std;

/*
求两个矩阵运算
*/

double*  MatrixMutl(double *a, double *b, int n);
double** MatrixMutl(double **a, double **b, int n);
double*  MatrixMutl(double **a, double *b, int n);
double*  MatrixMutl(double *a, int n, double t);
double*  MatrixAdd(double *a, double *b, int n);
double** MatrixAdd(double **a, double **b, int n);
double*  MatrixMinus(double *a,double *b,int n);
double** MatrixMinus(double **a, double **b, int n);
double** MatrixTranspose(double **a, int n);
double*  MatrixTranspose(double *a, int n);
/*
求两个向量的范数
*/
double Norm(double *a, int n);
//a-b的范数
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

/*
预处理共轭梯度法
*/


double *Gradient(int n,double *x0,double **a,double *b,double *c,int N,double TOL);
void testGradient();