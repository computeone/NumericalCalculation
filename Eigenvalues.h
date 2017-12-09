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

/*
Wielandt������
*/

void Wielandt(double *v, double **a, int n, double r0, double TOL, int N);
void testWielandt();

/*
Householder�任�㷨
*/
void Householder(double **a, int n);
void testHouseholder();


/*
QR�ֽ��㷨
*/

double QR(double **a, double *b, int n,double TOL,int M);
void testQR();
