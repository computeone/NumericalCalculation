#include "Integration.h"
#include <iostream>
#include <iomanip>
using namespace std;

double CompoundSimpson(double a, double b, int n,double(*function)(double)){
	double h = (b - a) / n;
	double XI0 = function(a) + function(b);
	//f(x2i-1)的和
	double XI1 = 0;
	//f(x2i)的和
	double XI2 = 0;
	for (int i = 1; i <n; i++){
		double X = a + i*h;
		if (i % 2 == 0){
			XI2 = XI2 + function(X);
			cout << "X2:" << XI2 << endl;
		}
		else{
			XI1 = XI1 + function(X);
			cout << "XI1:" << XI1 << endl;
		}
	}
	double XI = h*(XI0 + 2 * XI2 + 4 * XI1) / 3;
	return XI;
}
double functionCompoundSimpson(double d){
	return exp(d);
}
void testCompoundSimpson(){
	double result = CompoundSimpson(0, 4, 10,functionCompoundSimpson);
	cout << "积分值为:" << result << endl;
}
double** Romberg(double a, double b, int n,double(*function)(double)){
	double** R = new double *[n];
	for (int i = 0; i < n; i++){
		R[i] = new double[n];
	}
	double h = b - a;
	R[0][0] = (h / 2)*(function(a) + function(b));
	
	for (int i = 1; i < n; i++){
		double sum = 0;
		for (int k = 1; k <=pow(2, i - 1); k++){
			sum = sum + function(a + (k - 0.5)*h);
		}
		//梯形方法近似
		R[i][0] = 0.5*(R[i-1][0] + h*sum);
		for (int j = 1; j <=i; j++){
			//外推
			R[i][j] = R[i][j-1] + (R[i][j-1] - R[i-1][j-1]) / (pow(4, j) - 1);
		}
		h = h / 2;
		
	}
	return R;
}

void testRomberg(){
	double **p=Romberg(0, 3.1415926, 6,functionRomberg);
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < i + 1; j++){
			cout << "R[" << i << "][" << j << "]:" << setprecision(8)<<p[i][j] << " ";
		}
		cout << endl;
	}
}
double functionRomberg(double x){
	return sin(x);
}

double AdaptiveIntegration(double a, double b, double tol, int N, double(*function)(double)){
	double app = 0;
	int i = 1;
	double toli = 10 * tol;
	double ai = a;
	double hi = (b - a) / 2;
	double FAI = function(a);
	double FCI = function(a + hi);
	double FBI = function(b);

	//对整个区间的Simpson方法近似
	double SI = hi*(FAI + 4 * FCI + FBI) / 3;
	double LI = 1;
	
	while (i > 0){
		double FD = function(ai + hi / 2);
		double FE = function(ai + 3 * hi / 2);
		//对子区间一半的Simpson方法近似
		double S1 = hi*(FAI + 4 * FD + FBI) / 6;
		double S2 = hi*(FCI + 4 * FE + FBI) / 6;

		double v1 = ai;
		double v2 = FAI;
		double v3 = FCI;
		double v4 = FBI;
		double v5 = hi;
		double v6 = toli;
		double v7 = SI;
		double v8 = LI;

		//删除层次
		i = i - 1;
		if (fabs(S1 + S2 - v7) < v6){
			cout << "达到精度" << endl;
			app = app + (S1 + S2);
		}
		else{
			if (v8 >= N){
				cout << "算法失败！" << endl;
				return -1;
			}
			else{
				cout << "继续细分" << endl;
				//增加一个层次
				//右半子区间的数据
				i = i + 1;
				ai = v1 + v5;
				FAI = v3;
				FCI = FE;
				FBI = v4;
				hi = v5 / 2;
				toli = v6 / 2;
				SI = S2;
				LI = v8 + 1;

				//左子区间的数据
				i = i + 1;
				ai = v1;
				FAI = v2;
				FCI = FD;
				FBI = v3;
				hi = hi/2;
				toli = toli/2;
				SI = S1;
				LI = LI;
			} 
		}
	}
	 
	return app;

}
void testAdaptiveIntegration(){
	cout << "积分值：" << AdaptiveIntegration(0, 3.1415926, 0.0001, 5, functionRomberg);
}

double DoubleIntegration(double a, double b, double c, double d, int m, int n,
	double(*function)(double, double)){
	
	double h = (b - a) / n;
	double J1 = 0; //首选项
	double J2 = 0;
	double J3 = 0;

	for (int i = 0; i < n+1; i++){
		double x = a + i*h;
		double HX = (d - c) / m;
		double K1 = function(x, c) + function(x, d);
		double K2 = 0;
		double K3 = 0;

		for (int j = 1; j < m; j++){
			double y = c + j*HX;
			double Q = function(x, y);
			if (j % 2 == 0){
				K2 = K2 + Q;
			}
			else{
				K3 = K3 + Q;
			}
		}
		double L = (K1 + 2 * K2 + 4 * K3)*HX / 3;
		if (i == 0 || i == n){
			J1 = J1 + L;
		}
		else if (i % 2 == 0){
			J2 = J2 + L;
		}
		else{
			J3 = J3 + L;
		}
		
	}
	double J = h*(J1 + 2 * J2 + 4 * J3) / 3;
	return J;
}
double functionDoubleIntegration(double x, double y){
	return log(x + 2 * y);
}

void testDoubleIntegration(){
	double result = DoubleIntegration(1.4, 2.0, 1.0, 1.5, 4, 4, functionDoubleIntegration);
	cout << "二重积分值：" << setprecision(16)<<result << endl;
}

//m为4,n为4
double GaussDoubleIntegration(double a, double b, double c, double d,int m,int n,
	double(*function)(double,double)){
	double h1 = (b - a) / 2;
	double h2 = (b + a) / 2;
	double J = 0;

	for (int i = 0; i < m ; i++){
		double JX = 0;
		double x = h1 * r4[i] + h2;
		double d1 = d;
		double c1 = c;
		double k1 = (d1 - c1) / 2;
		double k2 = (d1 + c1) / 2;

		for (int j = 0; j < n; j++){
			double y = k1 * r4[j] + k2;
			double Q = function(x, y);
			JX = JX + c4[j] * Q;
		}
		J = J + c4[i] * k1*JX;
	}
	J = h1*J;
	return J;
}
void testGaussDoubleIntegration(){
	double result = GaussDoubleIntegration(1.4, 2.0, 1.0, 1.5, 4,4,functionDoubleIntegration);
	cout << "高斯二重积分求值：" << setprecision(16)<<result << endl;
}
double GaussTrebleIntegration(double a, double b, double c, double d, double e, double f,
	int m, int n, int p, double(*function)(double, double, double)){

	double h1 = (b - a) / 2;
	double h2 = (b + a) / 2;
	double J = 0;
	for (int i = 0; i < m ; i++){
		double JX = 0;
		double x = h1 * r4[i] + h2;
		double d1 = d;
		double c1 = c;
		double k1 = (d1 - c1) / 2;
		double k2 = (d1 + c1) / 2;

		for (int j = 0; j < n ; j++){
			double JY = 0;
			double y = k1 * r4[j] + k2;
			double e1 = e;
			double f1 = f;
			double l1 = (f1 - e1) / 2;
			double l2 = (f1 + e1) / 2;

			for (int k = 0; k < p; k++){
				double z = l1 * r4[k] + l2;
				double Q = function(x, y, z);
				JY = JY + c4[k] * Q;
			}
			JX = JX + c4[j] * JY;
		}
		J = J + c4[i] * k1*JX;
	}
	J = h1*J;
	return J;
}

