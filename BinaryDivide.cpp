#include <iostream>
#include <limits>
#include <iomanip>
#include "BinaryDivide.h"
using namespace std;


//f(x)=x^3+4*x^2-10  f(x)=0的根
double binaryDivide(double a,double b,double tol){
	cout << "二分法求根开始:" << endl;
	double result,p,fp;
	while (true){
		p = (a + b) / 2.0;
		fp = powf(p, 3) + 4 * powf(p, 2) - 10.0;

		cout <<setprecision(16)<<"p:" << p << endl;
		cout << setprecision(16)<<"fp:" << fp << endl<<endl;

		if (fp == 0){
			result = p;
			break;
		}
		if (fp < 0){
			a = p;
		}
		else{
			b = p;
		}
		
		//判断精度是否足够
		if (fabsf(fp) < tol){
			result = p;
			break;
		}

	}
	return result;
	cout << "二分法求根结束：" << endl;
}

double binaryDivide(double a, double b, double tol, double(*function)(double)){
	cout << "二分法求根开始：" << endl;

	double result, p, fp;
	while (true){
		p = (a + b) / 2.0;
		fp = function(p);

		cout << setprecision(16)<<"p:" << p << endl;
		cout << setprecision(16)<<"fp:" << fp << endl<<endl;
		if (fp == 0){
			result = p;
			break;
		}
		if (fp < 0){
			a = p;
		}
		else{
			b = p;
		}


		if (fabs(fp) < tol){
			result = p;
			break;
		}
	}

	return result;
	cout << "二分法求根结束:" << endl;
}