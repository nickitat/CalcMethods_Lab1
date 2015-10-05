#define _CRT_SECURE_NO_WARNINGS

#include "calculation.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;
using namespace Garbage;

double a;

matrix build_matrix(int N, int Nt, double h, double tau) {
	matrix A(N + 1, N + 1);
	for (int i = 0; i < N + 1; ++i) {
		if (i == 0 || i == N) {
			A[i][i] = 1;
		}
		else {
			A[i][i - 1] = -a*tau / (h*h);
			A[i][i] = 2 * a*tau / (h*h) + 1;
			A[i][i + 1] = -a*tau / (h*h);
		}
	}
	return A;
}

double f(double x, double t) {
	return x + 2 * t - exp(x) + a*(12 * x*x - 6 * x + t*exp(x));
}

double mu(double x) {
	return -x*x*x*x + x*x*x;
}

double u(double x, double t) {
	return -x*x*x*x + x*x*x + t*x + t*t - t*exp(x);
}

void calculate(int N, int Nt, double a) {
	double h = 1.0 / N;
	double tau = 1.0 / Nt;
	assert(2 * a*a*tau / (h*h) <= 1);
	cout << tau + h*h << endl;
	matrix A = build_matrix(N, Nt, h, tau);
	matrix Ainv = calculation::inverse_matrix::LU_method(A);
	matrix Uprev(N + 1, 1);
	for (int k = 0; k < N + 1; ++k) {
		Uprev[k][0] = mu(k*h);
	}
	long double triple_penetration = -1e20;
	for (int j = 1; j < Nt + 1; ++j) {
		matrix F(N + 1, 1);
		for (int k = 0; k < N + 1; ++k) {
			F[k][0] = tau*f(k*h, tau*(j));
		}
		for (int k = 0; k < N + 1; ++k) {
			Uprev[k][0] += F[k][0];
		}
		Uprev = Ainv * Uprev;
		Uprev[0][0] = u(0, tau*j);
		Uprev[N][0] = u(1, tau*j);
		long double penetration = -1e20;
		for (int k = 0; k < N + 1; ++k) {
			penetration = max(penetration, abs(Uprev[k][0] - u(h*k, tau*j)));
			//cout << Uprev[k][0] << " " << u(h*k, tau*j) << endl;
		}
		//cout << penetration << endl;
		//cout << endl << endl;
		triple_penetration = max(triple_penetration, penetration);
	}
	cout << endl << endl << triple_penetration << endl;
}

int main() {
	freopen("output.txt", "w", stdout);
	calculate(20, 20, 0.022);
	return 0;
}