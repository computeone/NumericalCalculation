#include "Horner.h"
#include <iostream>

using namespace std;


pair<double,double> Horner(double x0){
	double a[5];
	a[4] = 2;
	a[3] = 0;
	a[2] = -3;
	a[1] = 3;
	a[0] = -4;

	double y = a[4];
	double z = a[4];
	for (int j = 1; j < 4; j++){
		y = x0*y + a[4 - j];
		z = x0*z + y; 
	}
	y = x0*y + a[0];	
	return pair<double,double>(y,z);
}

pair<double,double> Horner(double x0,int n,double a[], double(*function)(double)){
	double y = a[n - 1];
	double z = a[n - 1];
	for (int j = 1; j < n - 1; j++){
		y = x0*y + a[n - 1 - j];
		z = x0*z + y;
	}
	y = x0*y + a[0];
	return pair<double,double>(y,z);
}

double testMuller(double x){
	return 16 * powf(x, 4) - 40 * powf(x, 3) + 5 * powf(x, 2) + 20 * x + 6;
}
double Muller(double x0, double x1, double x2, double tol, int n, double(*function)(double)){

	double h1 = x1 - x0;
	double h2 = x2 - x1;
	double p1 = (function(x1) - function(x0)) / h1;
	double p2 = (function(x2) - function(x1)) / h2;
	double d = (p2 - p1) / (h2 + h1);
	cout << d << endl;
	int i = 3;
	cout << "Muller算法开始：" << endl;
	while (i <= n){
		double b = p2 + h2*d;

		//可能复数运算
		double D = sqrtf(powf(b, 2) - 4 * function(x2)*d);
		double E;
		if (fabsf(b - D) < fabsf(b + D)){
			E = b + D;
		}
		else{
			E = b - D;
		}
		double h = -2 * function(x2) / E;
		double p = x2 + h;
		cout << "p:" << p << endl;
		if (fabs(h) < tol){
			cout << "Muller算法结束!" << endl;
			return p;
		}
		//设置下一次迭代
		x0 = x1;
		x1 = x2;
		x2 = p;
		h1 = x1 - x0;
		h2 = x2 - x1;
		p1 = (function(x1) - function(x0)) / h1;
		p2 = (function(x2) - function(x1)) / h2;
		d = (p2 - p1) / (h2 + h1);
		i = i + 1;
	}
	cout << "Muller算法失败！" << endl;
}