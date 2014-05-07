#include "ODE.h"
#include <iostream>
#include <iomanip>
using namespace std;


vector<pair<double,double>> Runge_Kutta(double a, double b, int N, double w, double(*function)(double, double)){

	vector<pair<double, double>> result;
	double h = (b - a) / N;
	double t = a;
	result.push_back(pair<double,double>(t, w));

	for (int i = 1; i < N+1; i++){
		double K1 = h*function(t, w);
		double K2 = h*function(t + h / 2, w + K1 / 2);
		double K3 = h*function(t + h / 2, w + K2 / 2);
		double K4 = h*function(t + h , w + K3);

		w = w + (K1 + 2 * K2 + 2 * K3 + K4)/6;
		t = a + i*h;
		 
		result.push_back(pair<double, double>(t, w));
	}

	return result;
}

void Runge_Kutta_Fehlberg(double a, double b, double TOL, double w, double hmax, double hmin,
	double(*function)(double, double)){
	double t = a;
	double h = hmax;
	bool FLAG = 1;
	cout << "t:" << setprecision(10) << t << " w:" << w << endl;
	
	while (FLAG == 1){
		double K1 = h*function(t, w);
		double K2 = h*function(t + h / 4.0, w + K1 / 4.0);
		double K3 = h*function(t + 3.0 / 8.0 * h, w + 3.0 / 32.0 * K1 + 9.0 / 32.0 * K2);
		double K4 = h*function(t + 12.0 / 13.0 * h, w + 1932.0 / 2197.0 * K1 - 7200.0 / 2197 * K2 + 7296.0 / 2197.0 * K3);
		double K5 = h*function(t + h, w + 439.0 / 216.0 * K1 - 8 * K2 + 3680.0 / 513.0 * K3 - 845.0 / 4104 * K4);
		double K6 = h*function(t + h / 2.0, w - 8.0 / 27.0 * K1 + 2.0 * K2 - 3544.0 / 2565.0 * K3 + 1859.0 / 4104.0 * K4 - 11.0 / 40.0 * K5);

		//R为误差
		double R = fabsf(K1 / 360.0 - 128.0 / 4275.0 * K3 - 2197.0 / 75240.0 * K4 + K5 / 50.0 + 2.0 / 55.0 * K6) / h;
		
		if (R <= TOL){
			//接受近似
			t = t + h;
			w = w + 25.0 / 216.0 * K1 + 1408.0 / 2565.0 * K3 + 2197.0 / 4104.0 * K4 - K5 / 5.0;
			cout << "t:" << setprecision(10) << t << "  w:" << setprecision(10) << w
				<< "  h:" << h <<" R:"<< setprecision(10)<<R<<endl;
		}
		//重新计算步长
		double o = 0.84*sqrtf(sqrtf((TOL / R)));
			
		if (o <= 0.1){
			h = 0.1*h;
		}
		else if (o >= 4){
			h = 4 * h;
		}
		else{
			h = o*h;
		}
		
		//判断算法是否失败
		if (h >= hmax){
			h = hmax;
		}
		if (t >= b){
			FLAG = 0;
		}
		else if (t + h > b){
			h = b - t;
		}
		else if (h < hmin){
			FLAG = 0;
			cout << "算法失败！" << endl;
			return;
		}

	}
	cout << "算法成功!" << endl;
}
void Adams(double a, double b, int N,double w, double(*function)(double, double)){

	double h = (b - a) / N;
	double t = a;
	double *tt = new double[4];
	double *ww = new double[4];
	tt[0] = t;
	ww[0] = w;
	cout << "t0:" << setprecision(10) << t << "  w0:" << setprecision(10) << w << endl;
	/*
	使用Runge-Kutta方法计算初始值
	*/
	for (int i = 1; i < 4; i++){
		double K1 = h*function(tt[i-1], ww[i-1]);
		double K2 = h*function(tt[i-1] + h / 2, ww[i-1] + K1 / 2);
		double K3 = h*function(tt[i-1] + h / 2, ww[i-1] + K2 / 2);
		double K4 = h*function(tt[i-1] + h, ww[i-1] + K3);

		ww[i] = ww[i-1] + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
		
		tt[i] = a + i*h;
		cout << "t:" << setprecision(10) << tt[i] << " w:"
			<< setprecision(10) << ww[i] << endl;
	}

	for (int i = 4; i < N + 1; i++){
		t = a + i*h;
		//预测wi
		w = ww[3] + h*(55.0*function(tt[3], ww[3]) - 59 * function(tt[2], ww[2]) +
			37 * function(tt[1], ww[1]) - 9 * function(tt[0], ww[0])) / 24.0;

		//校正wi
		w = ww[3] + h*(9 * function(t, w) + 19 * function(tt[3], ww[3]) - 5 * function(tt[2], ww[2])
			+ function(tt[1], ww[1])) / 24.0;

		cout << "t:" << setprecision(10) << t << " w:" << setprecision(10)
			<< w << endl;
		//准备下一次迭代
		for (int j = 0; j < 3; j++){
			tt[j] = tt[j + 1];
			ww[j] = ww[j + 1];
		}
		tt[3] = t;
		ww[3] = w;
	}

}

void Extrapolation(double a, double b, double w, double TOL, double hmax, double hmin,
	double(*function)(double, double)){
	double NK[] = { 2,4, 6, 8, 12, 16, 24, 32 };
	
	double TO = a;
	double WO = w;
	double h = hmax;
	bool FLAG = 1;
	double *y = new double[8];
	double **Q = new double*[8];
	for (int i = 0; i < 8; i++){
		Q[i] = new double[8];
	}

	for (int i = 1; i < 8; i++){
		for (int j = 1; j < i; j++){
			//Qij=hi^2/hi+1^2
			Q[i][j] = powf(NK[i + 1] / NK[j],2);
		}
	}

	while (FLAG == 1){
		int k = 1;
		bool NFLAG = 0;
		while (k <= 8 && NFLAG == 0){
			double HK = h / NK[k];
			double T = TO;
			double W2 = WO;
			double W3 = W2 + HK*function(T, W2);//Euler法的第一步
			T = TO + HK;

			for (int j = 1; j < NK[k]; j++){
				double W1 = W2;
				W2 = W3;
				W3 = W1 + 2 * HK*function(T, W3);//中点法
				T = TO + (j + 1)*HK;
			}

			//计算yk,1的端点校正
			double yk = (W3 + W2 + HK*function(T, W3)) / 2;
			//yk-1=yk-1,1,yk-2=yk-1,2....y1=yk-1,k-1
			if (k >= 2){
				int j = k;
				double v=y[k-1];
				while (j >= 2){
					/*
					外推法计算(yj-1=yk,k-j+2, yj-1=(hj-1^2*yj-hk^2*yj-1)/(hj-1^2-hk^2)
					*/
					y[j-1]= (y[j]-y[j-1])/(Q[k-1][j-1]-1);

					j = j - 1;  

				}
				if (fabsf(y[k-1]<v) < TOL){
					//y1被接受作为新的w
					NFLAG = 1;
				}
			}
			k = k + 1;		
		}

		k = k - 1;
		if (NFLAG == 0){
			//结果被拒绝
			h = h / 2;//w的新值被拒绝，减小h
			if (h < hmin){
				cout << "hmin exceeded!" << endl;
			}
			FLAG = 0;
		}
		//结果被接受
		else{
			WO = y[k-1];
			TO = TO + h;
			cout << "t:" << setprecision(10) << TO << "  w:" << setprecision(10) << WO
				<< " h:" << setprecision(10) << h << endl;
			if (TO >= b){
				FLAG = 0;
				cout << "算法成功完成！" << endl;
			}
			else if (TO + h > b){
				h = b - TO;
			}
			else if (k <= 3 && h < 0.5*hmax){
				h = 2 * h;
			}

		}
	}


}
void Equations_Runge_Kutta(double a, double b, int m, int N, double *w0, 
	double(**function)(double,double*)){

	double h = (b - a) / N;
	double t = a;
	
	double *k1 = new double[m];
	double *k2 = new double[m];
	double *k3 = new double[m];
	double *k4 = new double[m];
	double *w = new double[m];
	double *w1 = new double[m];
	double *w2 = new double[m];
	double *w3 = new double[m];

	for (int j = 0; j < m; j++){
		w[j] = w0[j];
	}
	cout << "t:" << t;
	for (int i = 0; i < m; i++){
		cout << "  w" << i << ":" << w0[i];
	}
	cout << endl;
	for (int i = 0; i < N; i++){

		for (int j = 0; j < m; j++){
			k1[j] = h*function[j](t, w);
		}

		for (int j = 0; j < m; j++){

			for (int p = 0; p < m; p++){
				w1[p] = w[p] + k1[p]/2;
			}

			k2[j] = h*function[j](t + h / 2,w1);
		}

		for (int j = 0; j < m; j++){
			for (int p = 0; p < m; p++){
				w2[p] = w[p] + k2[p] / 2;
			}
			k3[j] = h*function[j](t + h / 2, w2);
		}

		for (int j = 0; j < m; j++){

			for (int p = 0; p < m; p++){
				w3[p] = w[p] + k3[p];
			}

			k4[j] = h*function[j](t + h , w3);
		}

		for (int j = 0; j < m; j++){
			w[j] = w[j] + (k1[j] + 2*k2[j] + 2*k3[j]+k4[j])/6;
		}
		
		t = a + (i+1)*h;

		cout << "t:" << t;
		for (int i = 0; i < m; i++){
			cout << "  w" << i << ":" << w[i];
		}
		cout << endl;
	}
}
double Runge_Kutta_Function(double x, double y){
	return y - powf(x, 2) + 1;
}
void testRunge_Kutta(){
	vector<pair<double, double>> result = Runge_Kutta(0, 2, 10, 0.5, Runge_Kutta_Function);
	cout << "Runge_Kutta执行结果:" << endl;
	int i = 0;
	for (pair<double, double> p : result){
		cout << "迭代次数:" << ++i << "  t:" << setprecision(10)<<p.first << 
			"   w:" << setprecision(10)<<p.second << endl;
	}
}

void testRunge_Kutta_Fehlberg(){
	Runge_Kutta_Fehlberg(0, 2, 0.00001, 0.5, 0.25, 0.01, Runge_Kutta_Function);
}

void testAdams(){
	Adams(0, 2, 10, 0.5, Runge_Kutta_Function);
}

void testExtrapolation(){
	Extrapolation(0, 2, 0.5, 0.00000001, 0.25, 0.01, Runge_Kutta_Function);
}

double function1(double t,double *a){
	return a[1];
}
double function2(double t,double *a){
	double result=exp(2 * t)*sin(t) - 2 * a[0] + 2 * a[1];
	return result;
}
void testEquations_Runge_Kutta(){
	double(*function[2])(double, double*);
	function[0] = function1;
	function[1] = function2;
	double *w0 = new double[2];
	w0[0] = -0.4;
	w0[1] = -0.6;
	Equations_Runge_Kutta(0,1,2,10,w0,function);
}
