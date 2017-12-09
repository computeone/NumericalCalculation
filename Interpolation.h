#include <iostream>
#include <vector>
using namespace std;

/*
��ֵ�����ʽ�ƽ�
*/

double** Neville(double x,vector<double> xn, vector<double> fx);

/*
ţ�ٲ��̲�ֵ��
*/

double** divideNewton(vector<double> xn, vector<double> fx);
/*
Hermite����ʽ��ֵ
*/

double** Hermite(vector<double> xn, vector<double> fx,vector<double> fxx);

/*
��Ȼ����������ֵ
*/

double** CubicSpline(vector<double> xn, vector<double> fx);

//�����ݲ�������������ֵ
void testCubicSpline();