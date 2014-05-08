#include <iostream>
#include "LinearEquations.h"

using namespace std;

double* GaussianElimination(double **A, int n){
	double *x = new double[n-1];
	int p = 0;
	for (int i = 0; i < n-1; i++){
		p = n+1;
		for (int j = i; j < n-1; j++){
			if (A[j][i] != 0){
				if (j<p){
					p = j;
				}
			}
		}

		if (p==n+1){
			cout << "没有唯一解！" << endl;
			return NULL;
		}
		if (p != i){
			/*
			Ep<->Ei
			*/
			for (int m = 0; m < n; m++){
				double tmp = A[p][m];
				A[p][m] = A[i][m];
				A[i][m] = tmp;
			}
		}

		for (int k = i + 1; k < n-1; k++){
			double m = A[k][i] / A[i][i];
			/*
			Ej-Mji*Ei=Ej
			*/
			for (int q = 0; q < n; q++){
				A[k][q] = A[k][q] - m*A[i][q];
			}

		}
	}
	if (A[n - 2][n - 2] == 0){
		cout << "没有唯一解!" << endl;
		return NULL;
	}

	for (int i = 0; i < n-1; i++){
		for (int j = 0; j < n; j++){
			cout << A[i][j] << "  ";
		}
		cout << endl;
	}
	/*
	开始向后代换
	*/
	x[n - 2] = A[n - 2][n-1] / A[n - 2][n - 2];
	for (int i = n - 3; i >= 0; i--){
		double s = 0;
		for (int j = i + 1; j < n-1; j++){
			s = s + A[i][j] * x[j];
		}
		cout << "s"<<s << endl;
		x[i] = (A[i][n-1] -s)/ A[i][i];
	}

	return x;
}

void testGaussianElimination(){
	double **a = new double*[5];
	a[0] = new double[5];
	a[0][0] = 1;
	a[0][1] = -1;
	a[0][2] = 2;
	a[0][3] = -1;
	a[0][4] = -8;

	a[1] = new double[5];
	a[1][0] = 2;
	a[1][1] = -2;
	a[1][2] = 3;
	a[1][3] = -3;
	a[1][4] = -20;

	a[2] = new double[5];
	a[2][0] = 1;
	a[2][1] = 1;
	a[2][2] = 1;
	a[2][3] = 0;
	a[2][4] = -2;

	a[3] = new double[5];
	a[3][0] = 1;
	a[3][1] = -1;
	a[3][2] = 4;
	a[3][3] = 3;
	a[3][4] = 4;
	double *x = GaussianElimination(a,5);
	if (x == NULL){
	}
	else{
		cout << "方程组的解：" << endl;
		for (int i = 0; i < 4; i++){
			cout <<"x["<<i+1<<"]:"<< x[i] << "  ";
		}
		cout << endl;
	}
	
}