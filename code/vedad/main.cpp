#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "JacobiRotation.h"

using namespace std;

void exA(int);
void iterationDependency();

int main() {

//	iterationDependency();
	exA(20);

	return 0;
}


void exA(int n) {
	
	int iterations;
	JacobiRotation(n, &iterations);

}

void iterationDependency() {

	int iterations;
	
	fstream outFile;
	outFile.open("data/iterationsData.dat", ios::out);

	for (int n=50; n <= 300; n+=50) {

		JacobiRotation(n,&iterations);
		outFile << n << " " << iterations << endl;

	}

	outFile.close();

}


