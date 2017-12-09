#include <iostream>
#include <vector>
const double r2[] = { 0.5773502692, -0.5773502692 };
const double r3[] = { 0.7745966692, 0.00000000, -0.7745966692 };
const double r4[] = { 0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116 };
const double c2[] = {1.00000000,1.00000000};
const double c3[] = {0.5555555556,0.8888888889,0.5555555556};
const double c4[] = {0.3478548451,0.6521451549,0.6521451549,0.3478548451};
using namespace std;

/*
���ֵ���ֵ�㷨
*/
/*
����Simpson����
*/
double CompoundSimpson(double a, double b, int n);
void testCompoundSimpson();
double function(double d);

/*
Romberg���Ʒ�ʵ��
*/

double** Romberg(double a, double b, int n,double(*function)(double));
void testRomberg();
double functionRomberg(double x);

/*
����Ӧ�����
*/

double AdaptiveIntegration(double a, double b, double tol, int N, double(*function)(double));
void testAdaptiveIntegration();
double functionAdaptive(double x);

/*
���ػ�����ֵ
*/

//Simpson���ػ���
double DoubleIntegration(double a, double b,double c,double d,int m,int n,
	double(*function)(double,double));
double functionDoubleIntegration(double x, double y);
void testDoubleIntegration();

//Gauss���ػ���
double GaussDoubleIntegration(double a,double b,double c,double d,int m,int n,
	double(*function)(double,double));
void testGaussDoubleIntegration();

//Gauss���ػ���
double GaussTrebleIntegration(double a, double b,double c,double d,double e,double f,
	int m, int n, int p,double(*function)(double, double,double));
void testGaussTrebleIntegration();

