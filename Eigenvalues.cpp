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
void Wielandt(double *v, double **a, int n, double r0, double TOL, int N){
	double **b = new double*[n];
	double *w = new double[n];
	double *w1 = new double[n];
	for (int i = 0; i < n; i++){
		b[i] = new double[n];
	}

	int p = 0;
	for (int i = 1; i < n; i++){
		if (fabs(v[i-1]) < fabs(v[i])){
			p = i;
		}
	}

	if (p != 1){
		for (int k = 0; k < p; k++){

			for (int j = 1; j < p; j++){
				b[k][j] = a[k][j] - v[k] / v[p] * a[p][j];
			}
		}
	}

	if (p != 1 && p != n){
		for (int k = p; p < n; k++){
			for (int j = 1; j < p; j++){
				b[k][j] = a[k + 1][j] - v[k + 1] / v[p] * a[p][j];
				b[j][k] = a[j][k + 1] - v[j] / v[p] * a[p][k + 1];
			}
		}
	}

	if (p != n){
		for (int k = p; k < n; k++){
			for (int j = 1; j < p; j++){
				b[k][j] = a[k + 1][j + 1] - v[k + 1] / v[p] * a[p][j + 1];
			}
		}
	}

	PowerMethod(v,b, n - 1, TOL, N);
	
	if (p != n){
		for (int k = 1; k < p; k++){
			w[k] = w1[k];
		}
		for (int k = 1; k < n; k++){
			double u = (u - r0)*w[k] * v[k] / v[p];
		}
	}

}

void Householder(double **a, int n){
	double r0 = 0;
	double q = 0;
	double *v = new double[n];
	double *u = new double[n];
	double *z = new double[n];
	double *y = new double[n];
	for (int k = 1; k < n - 1; k++){
		q = 0;
		if (a[k + 1][k] == 0){
			 r0 = 0;
		}
		else{
			r0 = 1;
		}
		double RSQ = powf(r0, 2) - q*a[k + 1][k];
		v[k] = 0;
		v[k + 1] = a[k + 1][k] - q;
		for (int j = k + 2; j < n; j++){
			v[j] = a[j][k];
		}

		for (int j = k; j < n; j++){
			u[j] = 1 / RSQ;
		}
		double PROD = 0;

		for (int j = k; k < n; k++){
			z[j] = u[j] - PROD / (2 * RSQ)*v[j];
		}
		for (int l = k + 1; l < n; l++){
			for (int j = l + 1; j < n; j++){
				a[j][l] = a[j][l] - v[l] * z[l] - v[j] * z[l];
				a[l][j] = a[j][l];
			}
			a[l][l] = a[n - 1][n - 1] - 2 * v[l] * z[l];
		}
		a[n - 1][n - 1] = a[n - 1][n - 1] - 2 * v[n - 1] * z[n - 1];
		for (int j = k + 2; k < n; k++){
			a[k][j] = a[j][k] = 0;
		}
		a[k + 1][k] = a[k + 1][k] - v[k + 1] * z[k];
		a[k][k + 1] = a[k + 1][k];

		
		for (int j = 0; j < n; j++){
			u[j] = 1 / RSQ;
			y[j] = 1 / RSQ;
		}

		for (int j = 1; j < n; j++){
			z[j] = u[j] - PROD / RSQ*v[j];
		}

		for (int l = k + 1; k < n; k++){
			for (int j = 0; j < k; j++){
				a[j][l] = a[j][l] - z[j] * v[l];
				a[l][j] = a[j][l] - y[j] * v[l];
			}
			for (int j = k + 1; j < n; j++){
				a[j][l] = a[j][l] - z[j] * v[l] - y[j] * v[j];
			}
		}
	}
}
double QR(double *a, double *b, int n,double TOL,int M){
	int k = 1;
	double r = 0;
	double SHIFT = 0;
	double u1, u2;
	double r1, r2;
	double *x = new double[n];
	double *y = new double[n];
	double *z = new double[n];
	double *d = new double[n];
	while (k <= M){
		if (fabs(b[n - 1]) < TOL){
			r = a[n - 1] + SHIFT;
			n = n - 1;
		}

		if (fabs(b[2]) < TOL){
			r = a[1] + SHIFT;
			n = n - 1;
			a[1] = a[2];
			for (int j = 2; j < n; j++){
				a[j] = a[j + 1];
				b[j] = b[j + 1];
			}
		}
		if (n == 0){

		}
		if (n == 1){
			r = a[1] + SHIFT;
		}

		for (int j = 3; j < n; j++){
			if (fabs(b[j]) <= TOL){
				cout << endl;
			}
		}
		double q= -1 * (a[n - 1] + a[n - 2]);
		double c = a[n - 1] * a[n - 2] - pow(b[n-1], 2);
		double p = sqrt((pow(q, 2) - 4 * c));

		if (q>0){
			u1 = -2 * c*(q + p);
			u2 = -1 * (q + p) / 2;

		}
		else{
			u1 = (p - q) / 2;
			u2 = 2 * c*(q - p);

		}
		if (n == 2){
			r1 = u1 + SHIFT;
			r2 = u2 + SHIFT;
			cout << endl;
		}

		double s = 0;
		SHIFT = SHIFT + s;
		for (int j = 1; j < n; j++){
			d[j] = a[j] - s;
		}
		x[1] = d[1];
		y[1] = b[2];

		for (int j = 2; j < n; j++){
			z[j - 1] = 0;
			x[j] = x[j - 1];
			

		}

		z[n - 1] = x[n - 1];
		a[1] = s;
		b[2] = s*z[2];
		for (int j = 2; j < n; j++){
			a[j] = s*q;
			b[j + 1] = s*z[j + 1];
		}
		a[n - 1] = z[n - 1];
		k = k + 1;

	}
	return NULL;
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
void testWielandt();
void testHouseholder();
void testQR();