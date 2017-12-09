#include <iostream>
#include <vector>
using namespace std;


/*
Runge-Kutta����
*/

vector<pair<double,double>> Runge_Kutta(double a, double b, int N, double w, double(*function)(double, double));
void testRunge_Kutta();
double Runge_Kutta_Function(double x, double y);


/*
Runge-Kutta-Fehlberg�������㷨
*/

void Runge_Kutta_Fehlberg(double a, double b, double TOL, double w, double hmax, double hmin
	, double(*function)(double,double));
void testRunge_Kutta_Fehlberg();

/*
Adams 4��Ԥ��У���㷨
*/

void Adams(double a,double b,int N,double w,double(*function)(double,double));
void testAdams();

/*
���Ʒ�
*/

void Extrapolation(double a, double b, double w, double TOL, double hmax,
	double hmin, double(*function)(double, double));
void testExtrapolation();

/*
΢�ַ������Runge-Kutta�㷨
*/

void Equations_Runge_Kutta(double a, double b, int m, int N, double *w,
	double(**function)(double...));
void testEquations_Runge_Kutta();

/*
���Դ�з�
*/

void LinearShootMethod(double a, double b, double r0, double r1, int N,
	double(*function1)(double, double,double), double(*function2)(double, double, double));
void testLinearShootMethod();