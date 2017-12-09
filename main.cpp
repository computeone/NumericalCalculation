#include <iostream>
#include <algorithm>
#include "main.h"
using namespace std;

double test(double p){
	return powf(p, 3) + 4 * powf(p, 2) - 10.0;
};

void testBinaryDivide(){
	cout << "��Ϊ��"<<binaryDivide(1.0, 2.0,0.000002,&test) << endl;
}
void testFixedIterator(){
	cout << "��Ϊ:" << fixedIterator(0.1, 0.0000001, 10) << endl;
}
void testNewTonIterator(){
	cout << "��Ϊ:" << NewTonIterator(0.5, 3.1415926/4.0,0.000001,10,test) << endl;
}
void testSteffense(){
	cout << "��Ϊ��" << Steffensen(1.5, 0.00000001, 10) << endl;
}
void testHorner(){
	cout << "��Ϊ:" <<Horner(-2).second << endl;
	cout << "Muller��Ϊ:" << Muller(0.5, 1.0, 1.5, 0.00001, 10, testMuller) << endl;
}
void testdivideNewton(){
	cout << "n�����ϵ����" << endl;
	vector<double> xn;
	xn.push_back(1.0);
	xn.push_back(1.3);
	xn.push_back(1.6);
	xn.push_back(1.9);
	xn.push_back(2.2);
	vector<double> fx;
	fx.push_back(0.7651977);
	fx.push_back(0.6200860);
	fx.push_back(0.4554022);
	fx.push_back(0.2818186);
	fx.push_back(0.1103623);

	cout << divideNewton(xn,fx)<< endl;
}
void testNeville(){
	cout << "��Ϊ:" << endl;
	vector<double> xn;
	xn.push_back(1.0);
	xn.push_back(1.3);
	xn.push_back(1.6);
	xn.push_back(1.9);
	xn.push_back(2.2);
	vector<double> fx;
	fx.push_back(0.7651977);
	fx.push_back(0.6200860);
	fx.push_back(0.4554022);
	fx.push_back(0.2818186);
	fx.push_back(0.1103623);
	cout << Neville(1.5, xn, fx)[4][4] << endl;
}
void testHermite(){
	vector<double> xn;
	vector<double> fx;
	vector<double> fxx;
	xn.push_back(8.3);
	xn.push_back(8.6);
	fx.push_back(17.56492);
	fx.push_back(18.50515);
	fxx.push_back(3.116256);
	fxx.push_back(3.151762);
	double **q=Hermite(xn, fx, fxx);
	for (int i = 0; i < 2*(xn.size()-1) + 1; i++){
		cout << "Q[" << i << "][" << i << "]:" << q[i][i] << endl;
	}
}
int main(int argc,char** argv){
	//testFixedIterator();
	//testNewTonIterator();
	//testSteffense();
	//testHorner();
	//testNeville();
	//testdivideNewton();
	//testHermite();
	//testCubicSpline();
	//testRomberg();
	//testAdaptiveIntegration();
	//testDoubleIntegration();
	//testGaussDoubleIntegration();
	//testRunge_Kutta();
	testRunge_Kutta_Fehlberg();
	double a = 3.0 / 4 * 5;
	cout << a << endl;
}