#include "NewTonIterator.h"
#include <iostream>
using namespace std;

double NewTonIterator(double p0,double p1, double tol,int n){
	double q0 = cosf(p0)-p0;
	double q1 = cosf(p1)-p0;
	double p;
	int i = 0;
	cout << "牛顿迭代法开始：" << endl;
	while (i < (n + 1)){
		p = p1 - (q1*(p1 - p0)) / (q1 - q0);
		cout << "迭代次数:" << i <<endl;
		cout << "p:" <<p<< endl;
		
		if (fabs(p - p1) < tol){
			cout << "牛顿迭代法结束" << endl;
			return p;
		}
		p0 = p1;
		p1 = p;
		q0 = q1;
		q1 = cosf(p)-p1;
		i++;
	}
	cout << "牛顿迭代法失败！" << endl;
}

double NewTonIterator(double p0, double p1, double tol, int n, double(*function)(double)){
	double q0 = function(p0);
	double q1 = function(p1);
	double p;
	int i = 0;
	cout << "牛顿迭代法迭代开始：" << endl;
	while (i < (n + 1)){
		p = p1 - (q1*(p1 - p0)) / (q1 - q0);
		cout << "p:" <<p<< endl;
		if (fabs(p - p1) < tol){
			cout << "牛顿迭代法结束：" << endl;
			return p;
		}
		p0 = p1;
		p1 = p;
		q0 = q1;
		q1 = function(p);
		i++;
	}
	cout << "牛顿迭代法失败！" << endl;
}