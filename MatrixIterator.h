#include <iostream>

using namespace std;

/*
��������������
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
�����������ķ���
*/
double Norm(double *a, int n);
//a-b�ķ���
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

/*
Ԥ�������ݶȷ�
*/


double *Gradient(int n,double *x0,double **a,double *b,double *c,int N,double TOL);
void testGradient();