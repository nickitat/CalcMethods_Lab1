#ifndef __CALCULATION__H__
#define __CALCULATION__H__

#include "matrix_calculation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <utility>

using namespace std;

typedef double(*fptr_d)(double x);
typedef double(*fptr_dd)(double x, double t);
typedef double(*fptr_ddd)(double x, double t, double a);

namespace Garbage {

	class Calculation {

		typedef vector<double> vec;
		typedef vector<vector<double>> matrix;

	public:

		Calculation(double a, fptr_ddd f, fptr_dd u, fptr_d mu, int N, int Nt)
			: a(a), f(f), u(u), mu(mu), N(N), Nt(Nt)
		{}

		void calculate() {
			double h = 1.0 / N;
			double tau = 1.0 / Nt;

			// Courant condition
			assert(2 * a*a*tau / (h*h) <= 1);

			//cout << tau + h*h << endl;

			double a1 = -a*tau / (h*h);
			double a2 = 1 + 2 * a*tau / (h*h);
			double a3 = -a*tau / (h*h);
			matrix A(3, vec(N - 1, 0));
			for (int i = 0; i < N - 1; ++i) {
				if (i) A[0][i] = a1;
				A[1][i] = a2;
				if (i + 1 < N - 1) A[2][i] = a3;
			}
			vec U(N + 1);
			for (int i = 0; i <= N; ++i) {
				U[i] = u(h*i, 0);
			}
			double max_error = -1;
			for (int j = 1; j <= Nt; ++j) {
				vec F(N - 1);
				for (int i = 0; i < N - 1; ++i) {
					F[i] = tau*f(h*(i + 1), tau*(j + 1), a);
					if (i == 0) F[i] += a*tau / (h*h) * U[0];
					if (i == N - 2) F[i] += a*tau / (h*h) * U[N];
				}
				for (int i = 0; i < N - 1; ++i) {
					F[i] += U[i + 1];
				}
				solve_system(A, F, U);
				U[0] = u(0, tau*(j + 1));
				U[N] = u(1, tau*(j + 1));
				/*for (int i = 0; i <= N; ++i) {
					cout << U[i] << " " << u(h*i, tau*(j)) << endl;
					}
					cout << endl << endl;*/
				max_error = max(max_error, calc_error(U, tau*(j+1), h));
			}
			cout << max_error << endl;
		}

		void solve_system(matrix A, vec F, vec &U) {
			for (int i = 1; i < N-1; ++i) {
				A[1][i] -= (A[0][i] / A[1][i - 1]) * A[2][i - 1];
				F[i] -= (A[0][i] / A[1][i - 1]) * F[i - 1];
				A[0][i] = 0;
			}
			U[N - 1] = F[N - 2] / A[1][N - 2];
			for (int i = N - 2; i > 0; --i) {
				U[i] = (F[i-1] - A[2][i-1] * U[i + 1]) / A[1][i-1];
			}
			//for (int i = 0; i < N - 1; ++i) {
			//	/*for (int j = 0; j < 3; ++j) {
			//		cout << A[j][i] << " ";
			//	}*/
			//	cout << F[i];
			//	cout << endl;
			//}
			//cout << endl << endl;
		}

	private:

		const double a;
		const int N, Nt;
		const fptr_d mu;
		const fptr_dd u;
		const fptr_ddd f;

		double calc_error(vec &U, double t, double h) {
			double error = 0;
			for (int i = 0; i < U.size(); ++i) {
				error += abs(U[i] - u(h*i, t));
			}
			return error / (N + 1);
		}

	};

}

#endif/*__CALCULATION__H__*/