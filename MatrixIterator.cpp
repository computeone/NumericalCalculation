#include <iostream>
#include <iomanip>
#include "MatrixIterator.h"

using namespace std;

double Norm(double *a, double *b, int n){
	double s = 0;
	for (int i = 0; i < n; i++){
		s = s + pow(a[i] - b[i], 2);
	}

	return sqrt(s);
}
double* Jacobi(int n, double *x0, double **a, int N, double TOL){
	double *x = new double[n-1];
	int k = 1;
	std::cout<< "第0次迭代：" << endl;
	for (int i = 0; i < n - 1; i++){
		std::cout << setprecision(10) << setiosflags(ios::showpoint) << x0[i] << "  ";
	}
	std::cout << endl;
	while (k <= N){


		/*
		xi迭代
		*/
		for (int i = 0; i < n-1; i++){
			double s = 0;
			for (int j = 0; j < n-1; j++){
				if (j != i){
					s = s + a[i][j] * x0[j];
				}
			}
			x[i] = (-1 * s + a[i][n - 1])/a[i][i];
		}


		std::cout << "第" << k << "次迭代：" << endl;
		for (int i = 0; i < n - 1; i++){
			std::cout << setprecision(10) << setiosflags(ios::showpoint) << x[i] << "  ";
		}
		std::cout << endl;
		//std::cout << "norm:" << Norm(x, x0, n-1) << endl;
		if (Norm(x, x0,n-1) < TOL){
			return x;
		}
		k = k + 1;

		for (int i = 0; i < n-1; i++){
			x0[i] = x[i];
		}
	}
	std::cout << "最大迭代次数超出！" << endl;
	return NULL;
}

double * Gauss_Seidel(int n, double *x0, double **a, int N, double  TOL){
	double *x = new double[n - 1];
	int k = 1;

	std::cout << "第0次迭代：" << endl;
	for (int i = 0; i < n - 1; i++){
		std::cout << setprecision(10) << setiosflags(ios::showpoint) << x0[i] << "  ";
	}
	std::cout << endl;
	while (k <= N){
		/*
		xi迭代计算
		*/
		for (int i = 0; i < n - 1; i++){
			double s = 0;
			for (int j = 0; j < n - 1; j++){
				if (j < i){
					s = s + a[i][j] * x[j];
				}
				if (j>i){
					s = s + a[i][j] * x0[j];
				}
			}

			x[i] = (-1 * s + a[i][n - 1]) / a[i][i];
		}

		std::cout << "第" << k << "次迭代：" << endl;
		for (int i = 0; i < n - 1; i++){
			std::cout << setprecision(10) << setiosflags(ios::showpoint) << x[i] << "  ";
		}
		std::cout << endl;


		if (Norm(x, x0,n-1) < TOL){
			return x;
		}

		k = k + 1;
		for (int i = 0; i < n - 1;i++){
			x0[i] = x[i];
		}
	}
	cout << "迭代次数超过最大次数！" << endl;
}
double * SOR(double w,int n, double *x0, double** a, int N, double TOL){
	double *x = new double[n - 1];
	int k = 1;

	std::cout << "第0次迭代：" << endl;
	for (int i = 0; i < n - 1; i++){
		std::cout << setprecision(10) << setiosflags(ios::showpoint) << x0[i] << "  ";
	}
	std::cout << endl;

	while (k <= N){
		/*
		xi迭代计算
		*/

		for (int i = 0; i < n - 1; i++){
			double s = 0;
			for (int j = 0; j < n - 1; j++){
				if (j < i){
					s = s + a[i][j] * x[j];
				}
				if (j>i){
					s = s + a[i][j] * x0[j];
				}
			}
			x[i] = (1 - w)*x0[i] + w*(-1 * s + a[i][n - 1]) / a[i][i];
		}

		std::cout << "第" << k << "次迭代：" << endl;
		for (int i = 0; i < n - 1; i++){
			std::cout << setprecision(10) << setiosflags(ios::showpoint) << x[i] << "  ";
		}
		std::cout << endl;

		if (Norm(x, x0, n - 1) < TOL){
			return x;
		}

		k = k + 1;
		for (int i = 0; i < n - 1; i++){
			x0[i] = x[i];
		}
	}
	cout << "迭代次数超过最大迭代次数！" << endl;
}
void testJacobi(){
	double **a = new double*[5];
	for (int i = 0; i < 5; i++){
		a[i] = new double[5];
	}
	a[0][0] = 10;
	a[0][1] = -1;
	a[0][2] = 2;
	a[0][3] = 0;
	a[0][4] = 6;

	a[1][0] = -1;
	a[1][1] = 11;
	a[1][2] = -1;
	a[1][3] = 3;
	a[1][4] = 25;

	a[2][0] = 2;
	a[2][1] = -1;
	a[2][2] = 10;
	a[2][3] = -1;
	a[2][4] = -11;

	a[3][0] = 0;
	a[3][1] = 3;
	a[3][2] = -1;
	a[3][3] = 8;
	a[3][4] = 15;

	double *x0 = new double[4];
	x0[0]=0;
	x0[1]=0;
	x0[2]=0;
	x0[3]=0;

	double *x = Jacobi(5, x0, a, 10, 0.01);
}
void testGauss_Seidel(){
	double **a = new double*[5];
	for (int i = 0; i < 5; i++){
		a[i] = new double[5];
	}
	a[0][0] = 10;
	a[0][1] = -1;
	a[0][2] = 2;
	a[0][3] = 0;
	a[0][4] = 6;

	a[1][0] = -1;
	a[1][1] = 11;
	a[1][2] = -1;
	a[1][3] = 3;
	a[1][4] = 25;

	a[2][0] = 2;
	a[2][1] = -1;
	a[2][2] = 10;
	a[2][3] = -1;
	a[2][4] = -11;

	a[3][0] = 0;
	a[3][1] = 3;
	a[3][2] = -1;
	a[3][3] = 8;
	a[3][4] = 15;

	double *x0 = new double[4];
	x0[0] = 0;
	x0[1] = 0;
	x0[2] = 0;
	x0[3] = 0;

	double *x = Gauss_Seidel(5, x0, a, 10, 0.0001);
}
void testSOR(){
	double **a = new double*[5];
	for (int i = 0; i < 5; i++){
		a[i] = new double[5];
	}
	a[0][0] = 10;
	a[0][1] = -1;
	a[0][2] = 2;
	a[0][3] = 0;
	a[0][4] = 6;

	a[1][0] = -1;
	a[1][1] = 11;
	a[1][2] = -1;
	a[1][3] = 3;
	a[1][4] = 25;

	a[2][0] = 2;
	a[2][1] = -1;
	a[2][2] = 10;
	a[2][3] = -1;
	a[2][4] = -11;

	a[3][0] = 0;
	a[3][1] = 3;
	a[3][2] = -1;
	a[3][3] = 8;
	a[3][4] = 15;

	double *x0 = new double[4];
	x0[0] = 0;
	x0[1] = 0;
	x0[2] = 0;
	x0[3] = 0;

	double *x = SOR(1.25,5, x0, a, 10, 0.0001);
}