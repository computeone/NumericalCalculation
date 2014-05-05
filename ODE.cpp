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