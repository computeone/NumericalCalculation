#include <iostream>
#include <iomanip>

using namespace std;

/*
Pade�ƽ�����
*/

void Pade(int m, int n, double *a,double(*function)(double));
void testPade();

/*
Chebyshev����ƽ�
*/

void Chebyshev(int m, int n, double *a,double(*function)(double));
void testChebyshev();

double padeFunction(double x);