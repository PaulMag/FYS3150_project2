#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "JacobiRotation.h"

using namespace std;
using namespace arma;

double offDiagonal(mat, int*, int*, int);

void JacobiRotation(int n, int *iterations) {
	
	vec d(n);
	vec rho(n);
	vec V(n);
	
	double rho_min = 0.0;
	double rho_max = 4.5;

	// Step size
	double h = (rho_max - rho_min)/(n);

	double e = -1./(h*h);

	for (int i=0; i < n; i++) {

		rho(i) = rho_min + (i+1) * h;
		V(i) = rho(i) * rho(i);
		d(i) = 2./(h * h) + V(i);

	}

	mat A = zeros<mat>(n,n);

	for (int i=0; i<n; i++) {

		for (int j=0; j<n; j++) {
			
			if (i==j) { A(i,j) = d(i); }
			
			else if (abs(i-j) == 1) {
				
				A(i,j) = e;
			
			}

		}

	}

	double epsilon = 1e-8;

	int k;
	int l;

	*iterations = 0;
	
	while (offDiagonal(A, &k, &l, n) > epsilon) {
	
		double tau, t, c, s;
		double a_ll, a_kk, a_kl, a_ik, a_il;

		a_ll = A(l,l);
		a_kk = A(k,k);
		a_kl = A(k,l);

		// Set the elements with indices k, l to zero
		// in order to solve the other equations
		A(k,l) = 0.0;
		A(l,k) = 0.0;


		tau = (a_ll - a_kk)/(2 * a_kl);

		if (tau > 0) {

			t = 1./(tau + sqrt(1 + tau * tau));

		}

		else {
			
			t = 1./(tau - sqrt(1 + tau * tau));
		}

		c = 1./sqrt(1 + t * t);
		s = t * c;

		// Update the other elements in A
		A(k,k) = (a_kk * c * c) - (2 * a_kl * c * s) + (a_ll * s * s);	
		A(l,l) = (a_ll * c * c) + (2 * a_kl * c * s) + (a_kk * s * s);
		
		for (int i=0; i < n; i++) {

			if (i != k and i != l) {

				// For easier readability
				a_ik = A(i,k);
				a_il = A(i,l);

				A(i,k) = a_ik * c - a_il * s;
				A(k,i) = A(i,k); // Symmetric case

				A(i,l) = a_il * c + a_ik * s;
				A(l,i) = A(i,l); // Symmetric case
			}

		}

		(*iterations)++;
//		cout << "New matrix:" << endl;
//		cout << A << endl;


	}
	
	cout << "Iterations: " << *iterations << endl;
	vec sortEigenvals = sort(A.diag());	// Sort by increasing eigenvalues

	for (int i=0; i<3; i++) {

		cout << "Eigenvalue no. " << i+1 << ":	" << sortEigenvals(i) << endl; 

	}

	

		

}

double offDiagonal(mat A, int *k, int *l, int n) {

	double maxElement;

	for (int i=0; i < n; i++) {

		for (int j=i+1; j < n; j++) {

			if (abs(A(i,j)) > maxElement) {
				
				maxElement = abs(A(i,j));
				*k = i;
				*l = j;

			}
			
		}
	}

	return maxElement;
}


