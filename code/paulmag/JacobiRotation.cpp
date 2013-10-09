#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

void JacobiRotation(mat& A, int n, mat& R) {

    int k, l;

    double a_ik, a_il;
    double a_ll, a_kl, a_kk; //, a_lk
    double b_ik, b_il;
    double tau, t, s, c;
    double r_ik, r_il;

    double eps = 1e-8;
    double maxA = 10 * eps;
    int iterations = 0;

    // Matrix for eigenvectors, begins as the unit matrix:
    R = eye<mat>(n, n);

    // The algorithm:
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
        //a_lk = A(l,k);
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
                a_ik = A(i,k);
                a_il = A(i,l);
                b_ik = a_ik * c - a_il * s;
                b_il = a_il * c + a_ik * s;
                A(i,k) = b_ik;
                A(k,i) = b_ik;
                A(i,l) = b_il;
                A(l,i) = b_il;
            }
            r_ik = R(i,k);
            r_il = R(i,l);
            R(i,k) = c * r_ik - s * r_il;
            R(i,l) = c * r_il + s * r_ik;
        }

        A(k,k) = a_kk * c*c - 2*a_kl * c * s + a_ll * s*s;
        A(l,l) = a_ll * c*c + 2*a_kl * c * s + a_kk * s*s;
        A(k,l) = 0.0;
        A(l,k) = 0.0;

        iterations++;
    }

    cout << "Iterations: " << iterations << endl;
}
