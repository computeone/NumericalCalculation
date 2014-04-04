#include "BinaryDivide.h"
#include <iostream>

using namespace std;

double test(double p){
	return powf(p, 3) + 4 * powf(p, 2) - 10.0;
};

int main(){
	cout << "¸ùÎª£º"<<binaryDivide(1.0, 2.0,0.000001,&test) << endl;
}