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
Runge-Kutta-Fehlberg������
*/

void Runge_Kutta_Fehlberg(double a, double b, double TOL, double w, double hmax, double hmin
	, double(*function)(double,double));
void testRunge_Kutta_Fehlberg();