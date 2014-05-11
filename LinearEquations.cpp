#include <iostream>
#include <iomanip>
#include "LinearEquations.h"

using namespace std;

double* GaussianElimination(double **A, int n)  {
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

double* GaussianPivotElimination(double **a, int n){
	double *x = new double[n];
	int *NROW = new int[n];
	for (int i = 0; i < n; i++){
		NROW[i] = i;
	}
	int p;
	for (int i = 0; i < n-1; i++){
		p = i;
		for (int j = i; j < n; j++){
			if (fabs(a[j][i])>fabs(a[p][i])){
				p = j;
			}
		}

		if (a[p][i]==0){
			cout << "没有唯一解！" << endl;
			return NULL;
		}

		if (NROW[i] != NROW[p]){
			int NCOPY = NROW[i];
			NROW[i] = NROW[p];
			NROW[p] = NCOPY;
		}
		
		for (int j = i + 1; j < n; j++){
			double m = a[NROW[j]][i] / a[NROW[i]][i];
			
			for (int q = 0; q < n+1; q++){
				a[NROW[j]][q] = a[NROW[j]][q] - m*a[NROW[i]][q];
			}
		}
		if (a[NROW[n-1]][n-1] == 0){
			cout << "没有唯一解！" << endl;
		}

		x[n - 1] = a[NROW[n - 1]][n] / a[NROW[n - 1]][n - 1];

		for (int i = n - 2; i <=0; i++){
			double s = 0;	
			for (int j = i + 1; j < n; j++){
				s = s + a[NROW[i]][j] * x[j];
			}
			
			x[i] = (a[NROW[i]][n]-s)/a[NROW[i]][i];
		}
	}
	return x;
}
void LU(int n, double **a){
	double **l = new double*[n];
	double **u = new double*[n];

	for (int i = 0; i < n; i++){
		l[i] = new double[n];
		u[i] = new double[n];
		for (int j = 0; j < n; j++){
			l[i][j] = 0;
			u[i][j] = 0;
		}
	}
	/*
	设置lii为1
	*/
	for (int i = 0; i < n; i++){
		l[i][i] = 1;
	}

	u[0][0] = a[0][0] / l[0][0];
	if (u[0][0] == 0){
		cout << "无法分解！" << endl;
		return;
	}

	for (int j = 1; j < n; j++){
		u[0][j] = a[0][j] / l[0][0];//U的第一行
		l[j][0] = a[j][0] / u[0][0];//L的第一列
	}

	for (int i = 1; i < n - 1; i++){

		double s = 0;
		for (int k = 0; k < i; k++){
			s = s + l[i][k] * u[k][i];
		}
		
		u[i][i] = (a[i][i] - s) / l[i][i];
		
		if (u[i][i] == 0){
			cout << "无法分解！" << endl;
			return;
		}

		for (int j = i + 1; j < n; j++){
			s = 0;
			for (int k = 0; k < i; k++){
				s = s + l[i][k] * u[k][j];
			}

			u[i][j] = (a[i][j] - s) / l[i][i];//U的第i行
			s = 0;
			for (int k = 0; k < i; k++){
				s = s + l[j][k] * u[k][i];
			}

			l[j][i] = (a[j][i] - s) / u[i][i];//L的第i列
		}

	}

	double s = 0;
	for (int k = 0; k < n - 1; k++){
		s = s + l[n - 1][k] * u[k][n - 1];
	}
	/*
	如果lnnUnn=0，则A=LU，但A是奇异的
	*/
	u[n - 1][n - 1] = (a[n - 1][n - 1]-s)/l[n-1][n-1];

	//打印输出结果
	cout << "L矩阵：" << endl;
	for (int i = 0; i < n; i++){

		for (int j = 0; j < n; j++){
			cout << setprecision(5)<<setiosflags(ios::showpoint)<<l[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "U矩阵：" << endl;
	for (int i = 0; i < n; i++){

		for (int j = 0; j < n; j++){
			cout << u[i][j] << "  ";
		}
		cout << endl;
	}
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
void testLU(){
	double **a = new double*[4];
	for (int i = 0; i < 4; i++){
		a[i] = new double[4];
	}
	a[0][0] = 1;
	a[0][1] = 1;
	a[0][2] = 0;
	a[0][3] = 3;

	a[1][0] = 2;
	a[1][1] = 1;
	a[1][2] = -1;
	a[1][3] = 1;

	a[2][0] = 3;
	a[2][1] = -1;
	a[2][2] = -1;
	a[2][3] = 2;

	a[3][0] = -1;
	a[3][1] = 2;
	a[3][2] = 3;
	a[3][3] = -1;
	LU(4, a);

}
void testGaussianPivotElimination(){
	double **a = new double*[3];
	for (int i = 0; i < 3; i++){
		a[i] = new double[3];
	}

	a[0][0]=0.003000;
	a[0][1]=59.14;
	a[0][2]=59.17;

	a[1][0] = 5.291;
	a[1][1] = -6.130;
	a[1][2] = 46.78;

	double *x=GaussianPivotElimination(a, 2);
	if (x == NULL){
	}
	else{
		cout << "方程组的解：" << endl;
		for (int i = 0; i < 2; i++){
			cout << "x[" << i + 1 << "]:" << x[i] << "  ";
		}
		cout << endl;
	}
}