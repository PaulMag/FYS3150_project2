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

    int n     = 5;
    //int nStep = n + 1;

    double rhoMin = 0;
    double rhoMax = 10;
    double h = (rhoMax - rhoMin) / (n + 1);
    double e = - 1 / (h*h); // non-diagonal matrix element

    // Create matrix A:
    mat A = zeros<mat>(n, n);
    for (int i=0; i<n-1; i++) {
        A(i,i) = 2 / (h*h) + pow(rhoMin + (i+1)*h, 2);
        A(i+1,i) = e;
        A(i,i+1) = e;
    }
    A(n-1,n-1) = 2 / (h*h) + pow(rhoMin + n*h, 2);

    cout << A << endl;

    int k = 0;
    int l = 0;

    double a_ll, a_lk, a_kl, a_kk;
    double tau, t, s, c;
    double b_ik, b_il;

    double eps = 1e-8;
    double maxA = 10 * eps;
    int iterations = 0;

    while (abs(maxA) > eps) {

        maxA = 0;
        for (int i=0; i<n-1; i++) {
            if ( maxA < abs(A(i+1,i)) ) {
                maxA = A(i+1,i);
                k = i+1;
                l = i;
            }
            if ( maxA < abs(A(i,i+1)) ) {
                maxA = A(i,i+1);
                k = i;
                l = i+1;
            }
        }

        a_kk = A(k,k);
        a_kl = A(k,l);
        a_lk = A(l,k);
        a_ll = A(l,l);

        tau = (a_ll - a_kk) / (2 * a_kl);
        t = - tau + sqrt(1 + tau*tau);

        c = 1 / sqrt(1 + t*t); // cos(theta)
        s = t * c;             // sin(theta)

        for (int i=0; i<n; i++) {
            if (i!=k and i!=l) {
                b_ik = A(i,k) * c - A(i,l) * s;
                b_il = A(i,l) * c - A(i,k) * s;
                A(i,k) = b_ik;
                A(k,i) = b_ik;
                A(i,l) = b_il;
                A(l,i) = b_il;
            }
        }
        A(k,k) = a_kk * c*c - 2*a_kl * c * s + a_ll * s*s;
        A(l,l) = a_ll * c*c - 2*a_kl * c * s + a_kk * s*s;
        A(k,l) = 0;
        A(l,k) = 0;

        iterations++;
        cout << maxA << endl;
    }

    cout << A << endl;
    cout << "Iterations: " << iterations << endl;
}
