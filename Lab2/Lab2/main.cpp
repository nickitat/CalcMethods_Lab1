#define _CRT_SECURE_NO_WARNINGS

#include "calculation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;
using namespace Garbage;

typedef double(*fptr_d)(double x);
typedef double(*fptr_dd)(double x, double t);
typedef double(*fptr_ddd)(double x, double t, double a);

class Calculation {

public:

	Calculation(double a, fptr_ddd f, fptr_dd u, fptr_d mu, int N, int Nt)
		: a(a), f(f), u(u), mu(mu), N(N), Nt(Nt)
	{}

	void calculate() {
		double h = 1.0 / N;
		double tau = 1.0 / Nt;

		assert(2 * a*a*tau / (h*h) <= 1);
		
		//cout << tau + h*h << endl;
		
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
				F[k][0] = tau*f(k*h, tau*j, a);
			}
			for (int k = 0; k < N + 1; ++k) {
				Uprev[k][0] += F[k][0];
			}
			Uprev = Ainv * Uprev;
			Uprev[0][0] = u(0, tau*j);
			Uprev[N][0] = u(1, tau*j);
			long double penetration = -1e20;
			for (int k = 0; k < N + 1; ++k) {
				penetration = max(penetration, abs(Uprev[k][0] - u(k*h, tau*j)));
				//cout << Uprev[k][0] << " " << u(h*k, tau*j) << endl;
			}
			//cout << penetration << endl;
			//cout << endl << endl;
			triple_penetration = max(triple_penetration, penetration);
		}
		cout << endl << endl << triple_penetration << endl;
	}

private:

	const double a;
	const int N, Nt;
	const fptr_d mu;
	const fptr_dd u;
	const fptr_ddd f;

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

};

double f(double x, double t, double a) {
	return x + 2 * t - exp(x) + a*(12 * x*x - 6 * x + t*exp(x));
}

double mu(double x) {
	return -x*x*x*x + x*x*x;
}

double u(double x, double t) {
	return -x*x*x*x + x*x*x + t*x + t*t - t*exp(x);
}

int main() {
	freopen("output.txt", "w", stdout);
	double a = 0.022;
	int N = 20;
	int Nt = 20;
	Calculation calc(a, f, u, mu, N, Nt);
	calc.calculate();
	return 0;
}