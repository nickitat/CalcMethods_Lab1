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

			matrix A = build_matrix(N, Nt, h, tau);
			matrix Ainv = matrix_calculation::inverse_matrix::LU_method(A);

			matrix Uprev = get_initial_U(h);

			long double maximal_error = -1e20;

			for (int j = 1; j < Nt + 1; ++j) {

				matrix F = get_F(j, tau, h);

				Uprev = Uprev + F;

				Uprev = Ainv * Uprev;

				Uprev[0][0] = u(0, tau*j);
				Uprev[N][0] = u(1, tau*j);

				/*************************************************************************/
				long double local_error = -1e20;
				for (int k = 0; k < N + 1; ++k) {
					local_error = max(local_error, abs(Uprev[k][0] - u(k*h, tau*j)));
					//cout << Uprev[k][0] << " " << u(h*k, tau*j) << endl;
				}
				//cout << penetration << endl;
				//cout << endl << endl;
				maximal_error = max(maximal_error, local_error);
				/*************************************************************************/
			}
			cout << endl << endl << maximal_error << endl;
		}

	private:

		const double a;
		const int N, Nt;
		const fptr_d mu;
		const fptr_dd u;
		const fptr_ddd f;

		matrix build_matrix(int N, int Nt, double h, double tau) {
			matrix A(N + 1, N + 1);
			double a1 = -a*tau / (h*h);
			double a2 = 2 * a*tau / (h*h) + 1;
			for (int i = 0; i < N + 1; ++i) {
				if (i == 0 || i == N) {
					A[i][i] = 1;
				}
				else {
					A[i][i - 1] = a1;
					A[i][i] = a2;
					A[i][i + 1] = a1;
				}
			}
			return A;
		}

		matrix get_F(int j, double tau, double h) {
			matrix F(N + 1, 1);
			for (int k = 0; k < N + 1; ++k) {
				F[k][0] = tau*f(k*h, tau*j, a);
			}
			return F;
		}

		matrix get_initial_U(double h) {
			matrix Uinit(N + 1, 1);
			for (int k = 0; k < N + 1; ++k) {
				Uinit[k][0] = mu(k*h);
			}
			return Uinit;
		}

	};

}

#endif/*__CALCULATION__H__*/