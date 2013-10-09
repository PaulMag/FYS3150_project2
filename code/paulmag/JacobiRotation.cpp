#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

mat JacobiRotation(mat A, int n) {

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
        s = t * c;               // sin(theta)

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

        iterations++;
    }

    cout << "Iterations: " << iterations << endl;
    return A;
}
