#include <iostream>
#include "Eigenvalues.h"
#include "MatrixIterator.h"
#include "LinearEquations.h"
#include <iomanip>

using namespace std;

double PowerMethod(double *x0, double **a, int n,double TOL, int N){
	int k = 1;
	double *x = new double[n];
	double *y = new double[n];
	/*
	求范数
	*/
	double xp = x0[0];
	for (int i = 1; i < n; i++){
		if (fabs(xp) < fabs(x0[i])){
			xp = x0[i];
		}
	}

	for (int i = 0; i < n; i++){
		x[i] = x0[i] / xp;
	}

	while (k <= N){
		MatrixMutl(a, x ,y, n);
		/*
		求y的范数
		*/
		double yp = y[0];
		for (int i = 1; i < n; i++){
			if (fabs(yp) < fabs(y[i])){
				yp = y[i];
			}
		}
		double u = yp;
		if (yp == 0){
			cout << "第次" << k << "迭代特征向量和特征值:" << endl;
			for (int i = 0; i < n; i++){
				cout << "x[" << i+1 << "]:" << x[i] << "  ";
			}
			cout << u << endl;
			cout << "A 有0特征值，选择新的x重试。" << endl;
			return 0;
		}
		double ERR = x[0] - y[0] / yp;
		for (int i = 1; i < n; i++){
			if (fabs(ERR) < fabs((x[i]-y[i]/yp))){
				ERR = x[i]-y[i]/yp;
			}
		}

		for (int i = 0; i < n; i++){
			x[i] = y[i] / yp;
		}

		cout << "第次"<<k<<"迭代特征向量和特征值:" << endl;
		for (int i = 0; i < n; i++){
			cout << "x[" << i+1 << "]:" << x[i] << "  ";
		}
		cout << u << endl;
		/*
		判断精度要求
		*/
		if (ERR < TOL){
			return u;
		}
		k = k + 1;
	}
	return 0;
	cout << "超过迭代次数。" << endl;
}

double SymmetryPowerMethod(double *x0, double **a, int n, double TOL, int N){
	int k = 1;
	double *x = new double[n];
	double *y = new double[n];
	/*
	求范数xp
	*/
	double xp = 0;
	for (int i = 0; i < n; i++){
		xp = xp + pow(x0[i], 2);
	}
	xp = sqrt(xp);
	for (int i = 0; i < n; i++){
		x[i] = x0[i] / xp;
	}

	while (k <= N){
		MatrixMutl(a, x, y, n);

		double u = 0;
		for (int i = 0; i < n; i++){
			u = u + x[i] * y[i];
		}
		double yp = 0;
		for (int i = 0; i < n; i++){
			yp = yp + pow(y[i],2);
		}
		yp = sqrtf(yp);
		if (yp == 0){
			cout << "第次" << k << "迭代特征向量和特征值:" << endl;
			for (int i = 0; i < n; i++){
				cout << "x[" << i + 1 << "]:" << x[i] << "  ";
			}
			cout << u << endl;
			cout << "A 有0特征值，选择新的x重试。" << endl;
			return 0;
		}
		double ERR = 0;
		for (int i = 0; i < n; i++){
			ERR = ERR + pow(x[i] - y[i] / yp,2);
		}
		ERR = sqrt(ERR);

		for (int i = 0; i < n; i++){
			x[i] = y[i] / yp;
		}

		cout << "第次" << k << "迭代特征向量和特征值:" << endl;
		for (int i = 0; i < n; i++){
			cout << "x[" << i + 1 << "]:" << x[i] << "  ";
		}
		cout << u << endl;

		if (ERR < TOL){
			return u;
		}
		k = k + 1;
	}
	return 0;
	cout << "超过最大迭代次数！" << endl;

}

//a大小应传入为n+1
double AntiPowerMethod(double *x0, double **a, int n, double TOL, int N){
	double *t = new double[n];
	double **I = new double*[n];
	double *x = new double[n];
	double *y = new double[n];

	MatrixMutl(a, x0, t,n);
	

	double q = 0;
	double qq = 0;
	for (int i = 0; i < n; i++){
		q = q + pow(x0[i], 2);
	}

	q = MatrixMutl(x0, t, n)/ q;

	//初始化单位矩阵
	for (int i = 0; i < n; i++){
		I[i] = new double[n];
	}
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (i == j){
				I[i][j] = 1;
			}
			else{
				I[i][j] = 0;
			}
		}
	}
	delete[] t;

	int k = 1;
	double xp = x0[0];
	for (int i = 1; i < n; i++){
		if (fabs(xp) < fabs(x0[i])){
			xp = x0[i];
		}
	}
	for (int i = 0; i < n; i++){
		x[i] = x0[i] / xp;
	}

	while (k <= N){
		y=GaussianPivotElimination(MatrixMinus(a, MatrixMutl(I, n, q),n),x,n);
		if (y == NULL){
			cout << "q特征值：" << q << endl;
			return q;
		}

		double yp = y[0];
		for (int i = 0; i < n; i++){
			if (fabs(yp) < fabs(y[i])){
				yp = y[i];
			}
		}
		double u = yp;
		double ERR = x[0]-y[0]/yp;
		for (int i = 0; i < n; i++){
			if (fabs(ERR) < fabs((x[i] - y[i] / yp))){
				ERR = x[i] - y[i] / yp;
			}
		}

		for (int i = 0; i < n; i++){
			x[i] = y[i] / yp;
		}

		cout << "第次" << k << "迭代特征向量和特征值:" << endl;
		for (int i = 0; i < n; i++){
			cout << "x[" << i + 1 << "]:" << x[i] << "  ";
		}
		cout << 1/u+q << endl;
		if (ERR < TOL){
			u = 1 / u + q;
			return u;
		}
		k = k + 1;
	}
	return 0;
	cout << "超过最大迭代次数！" << endl;
	

}
void testPowerMethod(){
	double **a = new double*[3];
	for (int i = 0; i < 3; i++){
		a[i] = new double[3];
	}

	a[0][0]=-4;
	a[0][1]=14;
	a[0][2] =0;

	a[1][0] = -5;
	a[1][1] = 13;
	a[1][2] = 0;

	a[2][0] = -1;
	a[2][1] = 0;
	a[2][2] = 2;

	double *x0 = new double[3];
	x0[0] = 1;
	x0[1] = 1;
	x0[2] = 1;
	PowerMethod(x0, a, 3, 0.00001, 20);


}
void testSymmetryPowerMethod(){
	double **a = new double*[3];
	for (int i = 0; i < 3; i++){
		a[i] = new double[3];
	}

	a[0][0] = 4;
	a[0][1] = -1;
	a[0][2] = 1;

	a[1][0] = -1;
	a[1][1] = 3;
	a[1][2] = -2;

	a[2][0] = 1;
	a[2][1] = -2;
	a[2][2] = 3;

	double *x0 = new double[3];
	x0[0] = 1;
	x0[1] = 0;
	x0[2] = 0;

	SymmetryPowerMethod(x0, a, 3, 0.0001, 10);
}

void testAntiPowerMethod(){
	double **a = new double*[3];
	for (int i = 0; i < 3; i++){
		a[i] = new double[3];
	}

	a[0][0] = -4;
	a[0][1] = 14;
	a[0][2] = 0;

	a[1][0] = -5;
	a[1][1] = 13;
	a[1][2] = 0;

	a[2][0] = -1;
	a[2][1] = 0;
	a[2][2] = 2;

	double *x0 = new double[3];
	x0[0] = 1;
	x0[1] = 1;
	x0[2] = 1;
	AntiPowerMethod(x0, a, 3, 0.0001, 10);
}