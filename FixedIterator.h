#include <iostream>

using namespace std;


//不动点迭代法p=g(p)
double fixedIterator(double p0, double tol,int N);
double fixedIterator(double p0, double tol,int N, double(*function)(double p));