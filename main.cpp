#include "BinaryDivide.h"
#include "FixedIterator.h"
#include "NewTonIterator.h"
#include <iostream>
#include <algorithm>

using namespace std;

double test(double p){
	return powf(p, 3) + 4 * powf(p, 2) - 10.0;
};

void testBinaryDivide(){
	cout << "根为："<<binaryDivide(1.0, 2.0,0.000002,&test) << endl;
}
void testFixedIterator(){
	cout << "根为:" << fixedIterator(0.1, 0.0000001, 10) << endl;
}
void testNewTonIterator(){
	cout << "根为:" << NewTonIterator(0.5, 3.1415926/4.0,0.000001,10,test) << endl;
}
int main(){
	testFixedIterator();
	testNewTonIterator();

}