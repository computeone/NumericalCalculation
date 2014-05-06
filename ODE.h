#include <iostream>
#include <vector>
using namespace std;


/*
Runge-Kutta方法
*/

vector<pair<double,double>> Runge_Kutta(double a, double b, int N, double w, double(*function)(double, double));
void testRunge_Kutta();
double Runge_Kutta_Function(double x, double y);


/*
Runge-Kutta-Fehlberg误差控制算法
*/

void Runge_Kutta_Fehlberg(double a, double b, double TOL, double w, double hmax, double hmin
	, double(*function)(double,double));
void testRunge_Kutta_Fehlberg();

/*
Adams 4阶预测校正算法
*/

void Adams(double a,double b,int N,double w,double(*function)(double,double));
void testAdams();

/*
外推法
*/

void Extrapolation(double a, double b, double w, double TOL, double hmax,
	double hmin, double(*function)(double, double));
void testExtrapolation();