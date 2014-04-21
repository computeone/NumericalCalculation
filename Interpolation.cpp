#include "Interpolation.h"
#include <iostream>

using namespace std;


double** Neville(double x,vector<double> xn, vector<double> fx){
	int size=xn.size();
	double **q = (double**)new double*[size];
	for (int i = 0; i < size; i++){
		q[i] = new double[size];
	}

	//计算Qnn的初始值
	for (int i = 0; i < size; i++){
		q[i][0] = fx[i];
	}

	//计算其他的Qnn
	for (int i = 1; i < size; i++){
		for (int j = 1; j <=i; j++){
			q[i][j] = ((x - xn[i - j])*q[i][j - 1] - (x - xn[i])*q[i - 1][j - 1]) / (xn[i] - xn[i - j]);
			cout << "Q["<<i<<"]"<<"["<<j<<"]:" << q[i][j] << endl;
		}
	}
	return q;
}

double** divideNewton(vector<double> xn, vector<double> fx){
	int size = xn.size();
	double **f = (double**)new double*[size];
	for (int i = 0; i < size; i++){
		f[i] = new double[size];
	}

	//计算Fn0的值
	for (int i = 0; i < size; i++){
		f[i][0] = fx[i];
	}

	//计算差商
	for (int i = 1; i < xn.size(); i++){
		for (int j = 1; j < i + 1; j++){
			f[i][j] = (f[i][j - 1] - f[i - 1][j - 1]) / (xn[i] - xn[i-j]);
			cout << "F[" << i << "][" << j << "]:" << f[i][j] << endl;
		}
	}
	return f;
}