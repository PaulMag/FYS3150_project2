#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

#include "time.h"
#include "JacobiRotation.cpp"

using namespace std;
using namespace arma;

/* I am beautiful.
 * Vedad is not beautiful, because he laughed at me.
 */

/* Stuff to remember:
 *
 * mat A = zeros<mat>(n, n);
 * vec a(n);
 */

int main() {

    int n = 150; // nStep = n+1

    double rhoMin = 0;
    double rhoMax = 5.; // 4.5 is ideal according to vedadh
    double rho, V;
    double h = (rhoMax - rhoMin) / (n + 1);
    double e = - 1 / (h*h); // non-diagonal matrix element

    double omega_r = 0.01; // 0.01, 0.5, 1, 5

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
        A(i,i) = 2. / (h*h) + V;
    }

    JacobiRotation(A, n, S);
    //vec eigenVals = eig_sym(A); // Armadillos built-in method

    vec eigenVals       = A.diag();
    vec eigenValsSorted = sort(eigenVals);
    int* eigenValsIndices = new int[n];

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (eigenValsSorted(i) == eigenVals(j)) {
                eigenValsIndices[i] = j;
                break;
            }
        }
    }

    cout << "Eigenvalues: " << endl
         << eigenValsSorted(0) << endl
         << eigenValsSorted(1) << endl
         << eigenValsSorted(2) << endl;

    ofstream outfile;

    // Choose wich file to write to:
    outfile.open("data/prob_no_coulomb.dat");
    //outfile.open("data/prob_coulomb0.dat");
    //outfile.open("data/prob_coulomb1.dat");
    //outfile.open("data/prob_coulomb2.dat");
    //outfile.open("data/prob_coulomb3.dat");
    //outfile.open("data/prob_coulomb4.dat");

    outfile << n << "," << rhoMax << "," << omega_r << endl;

    for (int i=0; i<n; i++) {
        // Write the eigenvectors (not normalized):
        outfile << S.col(eigenValsIndices[0])(i) << ","
                << S.col(eigenValsIndices[1])(i) << ","
                << S.col(eigenValsIndices[2])(i) << endl;
    }

    outfile.close();
    delete [] eigenValsIndices;

    /* Results:
     * When nStep \approx 168 then 3 lowest eigenvalues are correct
     * with 4 leading digits, rhoMax = 4.5.
     */
}
