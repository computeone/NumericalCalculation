#include <iostream>

using namespace std;;


//Horner�㷨�������ʽp(x)=AnX^n+An-1X^n-1+....a1X+a0=(x-x0)Q(x)+b0��ֵ�͵���

pair<double,double> Horner(double x0);
pair<double,double> Horner(double x0, int n,double a[],double(*function)(double));

double testMuller(double x);
// Muller�㷨�������ʽ��m�ظ�
double Muller(double x0, double x1, double x2, double tol, int n, double(*function)(double));
