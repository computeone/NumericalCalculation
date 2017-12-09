#include <iostream>
#include <vector>

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