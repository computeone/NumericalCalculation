#include <iostream>
#include "Approximation.h"
#include "LinearEquations.h"
using namespace std;

void Pade(int m, int n,double *a,double(*function)(double)){
	int N = m + n;
	double q0 = 1;
	double p0 = a[0];
	double **b = new double*[N + 1];
	for (int i = 0; i < N+1; i++){
		b[i] = new double[N+1];
		
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j <=i - 1; j++){
			if (j <= n){
				b[i][j] = 0;
			}
		}
		if (i <= n){
			b[i][i] = 1;
		}

		for (int j = i + 1; j < N; j++){
			b[i][j] = 0;
		}
		for (int j = 0; j <=i; j++){
			if (j < m){
				b[i][n + j] = -1 * a[i - j];
			}
		}
		for (int j = n + i+1; j < N; j++){
			b[i][j] = 0;
		}
		b[i][N] = a[i+1];
	}
	/*
	使用部分选主元法求解线性方程组
	*/
	double *x = GaussianPivotElimination(b, N);
	if (x == NULL){
		cout << "没有求得系数！" << endl;
	}
	else{
		cout << "有理逼近系数：" << endl;
		for (int i = 0; i < N; i++){
			cout << "x[" << i + 1 << "]:" << x[i] << "  ";
		}
		cout << endl;
	}
}
double padeFunction(double x){
	return exp(-1 * x);
}
void testPade(){
	double *a = new double[3];
	a[0] = 1;
	a[1] = -1;
	a[2] = 0.5;
	a[3] = -1.0/6;
	a[4] = 1.0/24;
	a[5] = -1.0/120;
	Pade(2, 3, a, padeFunction);
}

void Chebyshev(int m, int n, double *a, double(*function)(double)){
	int N = m + n;
	double q0 = 1;
	double **b = new double*[N + 1];
	for (int i = 0; i < N + 1; i++){
		b[i] = new double[N + 1];
	}
	/*
	构建矩阵B
	*/
	for (int i = 0; i < N; i++){
		for (int j = 0; j < i; j++){
			if (j <= n){
				b[i][j] = 0;
			}
		}
		if (i <= n){
			b[i][i] = 0;
		}
		for (int j = i + 1; j < n; j++){
			b[i][j] = 0;
		}
		for (int j = n; j < N; j++){
			if (i != 0){
				b[i][j] = -1 * (a[i + j - n] + a[abs(i - j + n)])/2;
			}
			else{
				b[i][j] = -1 * a[j-n]/2;
			}
		}
		if (i != 0){
			b[i][N] = a[i];
		}
		else{
			b[i][N] = a[i] / 2;
		}
	}
	/*
	利用选主元法求解线性方程组
	*/

	double *x=GaussianPivotElimination(b, N);
	if (x == NULL){
	}
	else{
		cout << "有理逼近系数：" << endl;
		for (int i = 0; i < N; i++){
			cout << "x[" << i + 1 << "]:" << x[i] << "  ";
		}
		cout << endl;
	}
}

void testChebyshev(){
	cout << "test" << endl;
	double *a = new double[6];
	a[0] = 1.266066;
	a[1]=-1.130318;
	a[2]=0.271495;
	a[3]=-0.044337;
	a[4]=0.005474;
	a[5]=-0.000543;
	Chebyshev(2, 3, a, padeFunction);
}