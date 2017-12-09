#include "Steffensen.h"
#include <iostream>

double Steffensen(double p0, double tol, int n){
	int i = 1;
	double p1, p2,p;
	cout << "Steffense���������㷨��ʼ��" << endl;
	while (i <= n){
		p1 = sqrtf(10/(p0+4));
		p2 = sqrtf(10/(p1+4));
		p = p0 - powf(p1 - p0, 2) / (p2 - 2 * p1 + p0);
		cout << "p:" << p << endl;
		if (fabs(p - p0) < tol){
			cout << "Steffense���������㷨����:" << endl;
			return p;
		}
		i += 1;
		p0 = p;
	}
	cout << "Stenffense���������㷨ʧ�ܣ�" << endl;
	return NULL;
}

double Steffensen(double p0, double tol, int n, double(*function)(double)){
	int i = 1;
	double p1, p2, p;
	cout << "Steffense���������㷨��ʼ��" << endl;
	while (i <= n){
		p1 = function(p0);
		p2 = function(p1);
		p = p0 - powf(p1 - p0, 2) / (p2 - 2 * p1 + p0);
		cout << "p:" << p << endl;
		if (fabsf(p - p0) < tol){
			cout << "Steffense���������㷨����:" << endl;
			return p;
		}
		i = i + 1;
		p0 = p;
	}

	cout << "Stenffense���������㷨ʧ�ܣ�" << endl;
	return NULL;
}