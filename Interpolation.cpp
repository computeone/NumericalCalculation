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

double** Hermite(vector<double> xn, vector<double> fx,vector<double> fxx){
	int size = xn.size();
	double **q = new double *[2*size+1];
	double *zn = new double[2 * size + 1];
	for (int i = 0; i < 2*size+1; i++){
		q[i] = new double[2*size+1];
	}

	
	for (int i = 0; i < size; i++){
		zn[2*i] = xn[i];
		zn[2*i + 1] = xn[i];
		q[2 * i][0] = fx[i];
		q[2 * i + 1][0] = fx[i];
		q[2 * i + 1][1] = fxx[i];
		if (i != 0){
			q[2 * i][1] = (q[2 * i][0] - q[2 * i - 1][0]) / (zn[2 * i] - zn[2 * i - 1]);
		}
	}

	for (int i = 2; i < 2 * size + 1; i++){
		for (int j = 2; j <= i; j++){
			q[i][j] = (q[i][j-1] - q[i - 1][j-1]) / (zn[i] - zn[i - j]);
		}
	}

	return q;

}


double** CubicSpline(vector<double> xn, vector<double> fx){
	int size = xn.size();
	double **p = new double *[size];
	for (int i = 0; i < size; i++){
		p[i] = new double[4];
	}
	//复制fx到p[][0]
	for (int i = 0; i < fx.size(); i++){
		p[i][0] = fx[i];
	}
	double *h = new double[size];
	double *l = new double[size];
	double *u = new double[size];
	double *z = new double[size];
	for (int i = 0; i < size - 1; i++){
		h[i] = xn[i + 1] - xn[i];
	}

	for (int i = 1; i < size - 1; i++){
		fx[i] = 3 / h[i] * (fx[i + 1] - fx[i]) - 3 / h[i - 1] * (fx[i] - fx[i - 1]);
	}

	l[0] = 1;
	u[0] = 0;
	z[0] = 0;

	for (int i = 1; i < size - 1; i++){
		l[i] = 2 * (xn[i + 1] - xn[i - 1]) - h[i - 1] * u[i - 1];
		u[i] = h[i] / l[i];
		z[i] = (fx[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	l[size - 1] = 1;
	z[size - 1] = 0;
	p[size - 1][2] = 0;

	for (int j = size - 2; j >= 0; j--){
		p[j][2] = z[j] - u[j] * p[j+1][2];
		p[j][1] = (fx[j + 1] - fx[j]) / h[j] - h[j] * (p[j + 1][2] + 2 * p[j][2]) / 3;
		p[j][3] = (p[j + 1][2] - p[j][2]) / (3 * h[j]);
	}

	return p;

}

void testCubicSpline(){

	vector<double> xn;
	vector<double> fx;
	xn.push_back(0.9);
	fx.push_back(1.3);
	xn.push_back(1.3);
	fx.push_back(1.5);
	xn.push_back(1.9);
	fx.push_back(1.85);
	xn.push_back(2.1);
	fx.push_back(2.1);
	
	double **p=CubicSpline(xn, fx);
	for (int i = 0; i < xn.size()-1; i++){
		cout << "-------------------" << endl;
		cout << "a:" << p[i][0] << endl;
		cout << "b:" << p[i][1] << endl;
		cout << "c:" << p[i][2] << endl;
		cout << "d:" << p[i][3] << endl;
		cout << "----------------------" << endl;
	}

}