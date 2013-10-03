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

    int n = 5;
    mat A = zeros<mat>(n, n);
    for (int i=0; i<n-1; i++) {
        // example values:
        A(i,i) = 5;
        A(i+1,i) = 4;
        A(i,i+1) = 3;
    }
    A(4,4) = 5;

    cout << A << endl;

    int k = 0;
    int l = 0;

    double a_ll, a_lk, a_kl, a_kk;
    double tau, t, s, c;
    double b_ik, b_il;

    double eps = 1e-8;
    double maxA = 10 * eps;
    int iterations = 0;

    while (maxA > eps) {

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

        for (int i=0; i<n and i!=k and i!=l; i++) {
            b_ik = A(i,k) * c - A(i,l) * s;
            b_il = A(i,l) * c - A(i,k) * s;
            A(i,k) = b_ik;
            A(k,i) = b_ik;
            A(i,l) = b_il;
            A(l,i) = b_il;
        }
        A(k,k) = a_kk * c*c - 2*a_kl * c * s + a_ll * s*s;
        A(l,l) = a_ll * c*c - 2*a_kl * c * s + a_kk * s*s;
        A(k,l) = 0;
        A(l,k) = 0;

        iterations++;
    }

    cout << A << endl;
    cout << "Iterations: " << iterations << endl;
}
