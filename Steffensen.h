#include <iostream>

using namespace std;


//Steffensen方法是具有二次收敛速度的算法,不动点迭代算法改进性
double Steffensen(double p0, double tol, int N);
double Steffensen(double p0, double tol, int N, double(*function)(double));