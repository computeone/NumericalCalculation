#include <iostream>

using namespace std;


//Steffensen�����Ǿ��ж��������ٶȵ��㷨,����������㷨�Ľ���
double Steffensen(double p0, double tol, int N);
double Steffensen(double p0, double tol, int N, double(*function)(double));