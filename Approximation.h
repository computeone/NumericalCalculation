#include <iostream>
#include <iomanip>

using namespace std;

/*
Pade逼近技术
*/

void Pade(int m, int n, double *a,double(*function)(double));
void testPade();

/*
Chebyshev有理逼近
*/

void Chebyshev(int m, int n, double *a,double(*function)(double));
void testChebyshev();

double padeFunction(double x);