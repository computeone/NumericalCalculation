#include "FixedIterator.h"
#include <iostream>
#include <iomanip>

using namespace std;


double fixedIterator(double p0, double tol,int N){

	int i = 0;
	double p;
	cout << "不动点迭代开始:" << endl;
	while (i < (N+1)){
		p = p0 - (powf(p0, 3) + 4 * powf(p0, 2) - 10.0) / (3 * powf(p0, 2) + 8 * p0);
		cout << setprecision(16)<<"p:" << p << endl << endl;
		if (fabs(p - p0) < tol){
			cout << "不动点迭代结束：" << endl;
			return p;
		}
		p0 = p;
	}
	cout << "不动点迭代结束:迭代失败!" << endl;
	return NULL;
}

double fixedIterator(double p0, double tol,int N, double(*function)(double p)){
	int i = 0;
	double p;
	cout << "不动点迭代开始:" << endl;
	while (i < (N + 1)){
		p = function(p0);
		cout << setprecision(16)<<"p:" << endl << endl;
		if (fabs(p - p0) < tol){
			cout << "不动点迭代结束：" << endl;
		}
		p0 = p;
	}
	cout << "不动点迭代结束：迭代失败！" << endl;
	return NULL;
}