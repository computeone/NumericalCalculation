#include <iostream>

using namespace std;


//�����������p=g(p)
double fixedIterator(double p0, double tol,int N);
double fixedIterator(double p0, double tol,int N, double(*function)(double p));