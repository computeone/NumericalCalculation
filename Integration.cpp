#include "Integration.h"
#include <iostream>

using namespace std;

double CompoundSimpson(double a, double b, int n,double(*function)(double)){
	double h = (b - a) / n;
	double XI0 = function(a) + function(b);
	//f(x2i-1)的和
	double XI1 = 0;
	//f(x2i)的和
	double XI2 = 0;
	for (int i = 1; i <n; i++){
		double X = a + i*h;
		if (i % 2 == 0){
			XI2 = XI2 + function(X);
			cout << "X2:" << XI2 << endl;
		}
		else{
			XI1 = XI1 + function(X);
			cout << "XI1:" << XI1 << endl;
		}
	}
	double XI = h*(XI0 + 2 * XI2 + 4 * XI1) / 3;
	return XI;
}

void testCompoundSimpson(){
	double result = CompoundSimpson(0, 4, 10,function);
	cout << "积分值为:" << result << endl;
}
double function(double d){
	return exp(d);
}