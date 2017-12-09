#include "NewTonIterator.h"
#include <iostream>
using namespace std;

double NewTonIterator(double p0,double p1, double tol,int n){
	double q0 = cosf(p0)-p0;
	double q1 = cosf(p1)-p0;
	double p;
	int i = 0;
	cout << "ţ�ٵ�������ʼ��" << endl;
	while (i < (n + 1)){
		p = p1 - (q1*(p1 - p0)) / (q1 - q0);
		cout << "��������:" << i <<endl;
		cout << "p:" <<p<< endl;
		
		if (fabs(p - p1) < tol){
			cout << "ţ�ٵ���������" << endl;
			return p;
		}
		p0 = p1;
		p1 = p;
		q0 = q1;
		q1 = cosf(p)-p1;
		i++;
	}
	cout << "ţ�ٵ�����ʧ�ܣ�" << endl;
}

double NewTonIterator(double p0, double p1, double tol, int n, double(*function)(double)){
	double q0 = function(p0);
	double q1 = function(p1);
	double p;
	int i = 0;
	cout << "ţ�ٵ�����������ʼ��" << endl;
	while (i < (n + 1)){
		p = p1 - (q1*(p1 - p0)) / (q1 - q0);
		cout << "p:" <<p<< endl;
		if (fabs(p - p1) < tol){
			cout << "ţ�ٵ�����������" << endl;
			return p;
		}
		p0 = p1;
		p1 = p;
		q0 = q1;
		q1 = function(p);
		i++;
	}
	cout << "ţ�ٵ�����ʧ�ܣ�" << endl;
}