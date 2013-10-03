#include <iostream>
#include <cmath>
#include <armadillo>
#include "JacobiRotation.h"

using namespace std;

void exA(int);

int main() {

	exA(100);

	return 0;
}


void exA(int n) {

	JacobiRotation(n);
}
