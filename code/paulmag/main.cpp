#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

/*
 *  I am beautiful.
 *  Vedad is not beautiful, because he laughed at me.
 */

/* Stuff to remember:
 *
 * mat A = zeros<mat>(n, n);
 * vec a(n);
 */

int main() {

    int n = 168;
    //int nStep = n + 1;

    double rhoMin = 0;
    double rhoMax = 4.5; // 4.5 is ideal according to vedadh
    double rho, V;
    double h = (rhoMax - rhoMin) / (n + 1);
    double e = - 1 / (h*h); // non-diagonal matrix element

    // Create matrix A:
    mat A = zeros<mat>(n, n);
    for (int i=0; i<n-1; i++) {
        rho = rhoMin + (i+1) * h;
        V = rho*rho;
        A(i,i) = 2 / (h*h) + V;
        A(i+1,i) = e;
        A(i,i+1) = e;
    }
    rho = rhoMin + n * h;
    V = rho*rho;
    A(n-1,n-1) = 2 / (h*h) + V;
    //A(n-1,n-1) = 2 / (h*h) + pow(rhoMin + n*h, 2);

    //cout << A << endl;

    int k, l;

    double a_ll, a_lk, a_kl, a_kk;
    double tau, t, s, c;
    double b_ik, b_il;

    double eps = 1e-8;
    double maxA = 10 * eps;
    int iterations = 0;

    while (maxA > eps) {

        maxA = 0;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                if (abs(A(i,j)) > maxA ) {
                    // Find max value except on diagonal.
                    maxA = abs(A(i,j));
                    k = i;
                    l = j;
                }
            }
        }

        a_kk = A(k,k);
        a_kl = A(k,l);
        a_lk = A(l,k);
        a_ll = A(l,l);

        tau = (a_ll - a_kk) / (2 * a_kl);
        if (tau > 0) {
            t = 1. / (tau + sqrt(1 + tau*tau));
        }
        else {
            t = 1. / (tau - sqrt(1 + tau*tau));
        }

        c = 1.0 / sqrt(1 + t*t); // cos(theta)
        s = t * c;             // sin(theta)

        for (int i=0; i<n; i++) {
            if (i!=k and i!=l) {
                b_ik = A(i,k) * c - A(i,l) * s;
                b_il = A(i,l) * c + A(i,k) * s;
                A(i,k) = b_ik;
                A(k,i) = b_ik;
                A(i,l) = b_il;
                A(l,i) = b_il;
            }
        }

        A(k,k) = a_kk * c*c - 2*a_kl * c * s + a_ll * s*s;
        A(l,l) = a_ll * c*c + 2*a_kl * c * s + a_kk * s*s;
        A(k,l) = 0.0;
        A(l,k) = 0.0;
        //cout << (a_kk - a_ll) * c*s + a_kl * (c*c - s*s) << endl;

        //cout << A << endl;
        iterations++;
        //cout << maxA << endl;
        //break;
    }

    //vec eigenVals = sort(eig_sym(A)); // Armadillos method

    cout << "Iterations: " << iterations << endl;
    cout << "Eigenvalues: " << endl;
    vec eigenVals = sort(A.diag());
    cout << eigenVals(0) << endl <<
            eigenVals(1) << endl <<
            eigenVals(2) << endl;

    /* When nStep \approx 168 then 3 lowest eigenvalues are correct
     * with 4 leading digits, rhoMax = 4.5.
     */

}
