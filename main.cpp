#include "BinaryDivide.h"
#include "FixedIterator.h"
#include "NewTonIterator.h"
#include "Steffensen.h"
#include "Horner.h"
#include <iostream>
#include <algorithm>

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

int main(){
	testFixedIterator();
	testNewTonIterator();
	testSteffense();
	testHorner();
}