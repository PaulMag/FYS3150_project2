#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

#include "time.h"
#include "JacobiRotation.cpp"

using namespace std;
using namespace arma;

/*
 *  I am beautiful.
 *  Vedad is not beautiful, because he laughed at me.
 */

/*
 * Stuff to remember:
 *
 * mat A = zeros<mat>(n, n);
 * vec a(n);
 */

int main() {

    int n = 200; // nStep = n+1

    double rhoMin = 0;
    double rhoMax = 4.5; // 4.5 is ideal according to vedadh
    double rho, V;
    double h = (rhoMax - rhoMin) / (n + 1);
    double e = - 1 / (h*h); // non-diagonal matrix element

    double omega_r = 0.0; // 0.01, 0.5, 1, 5

    // Create matrix A:
    mat A = zeros<mat>(n, n);
    mat S(n,n); // for keeping eigenvectors

    for (int i=0; i<n-1; i++) {
        A(i+1,i) = e;
        A(i,i+1) = e;
    }

    for (int i=0; i<n; i++) {
        rho = rhoMin + (i+1) * h;
        // Choose wich potential to use:
        V = rho*rho;
        //V = omega_r*omega_r * rho*rho + 1./rho;
        A(i,i) = 2 / (h*h) + V;
    }

    JacobiRotation(A, n, S);
    //vec eigenVals = sort(eig_sym(A)); // Armadillos method

    cout << "Eigenvalues: " << endl;
    vec eigenVals = sort(A.diag());
    cout << eigenVals(0) << endl
         << eigenVals(1) << endl
         << eigenVals(2) << endl;

    ofstream outfile;
    outfile.open ("data/prob_no_coulomb_0.dat");
    for (int i=0; i<n; i++) {
        outfile << pow(S.col(0)(i), 2) << endl;
    }
    outfile.close();

    /*
     * When nStep \approx 168 then 3 lowest eigenvalues are correct
     * with 4 leading digits, rhoMax = 4.5.
     */

}
