#ifndef __CALCULATION__H__
#define __CALCULATION__H__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

typedef double(*fptr_d)(double x);
typedef double(*fptr_dd)(double x, double t);
typedef double(*fptr_ddd)(double x, double t, double a);

namespace Garbage {

	class Calculation {

		typedef vector<double> vec;
		typedef vector<vector<double>> matrix;

	public:

		Calculation(fptr_d phi, fptr_dd C, fptr_dd u, int N, int Nt)
			: phi(phi), C(C), u(u), N(N), Nt(Nt)
		{}

		void calculate() {
			double h = 1.0 / N;
			double tau = 1.0 / Nt;

			vector<double> Uprev = getInitU(h);
			double maxError = -1;
			for (int j = 1; j <= Nt; ++j) {
				vector<double> Unew = solveSystem(Uprev, tau, h, j);
				double error = 0;
				for (int i = 0; i <= N; ++i) {
					error += abs(Unew[i] - u(h*i, tau*j));
				}
				maxError = max(maxError, error / (N + 1));
				Uprev = Unew;
			}
			cout << maxError << endl;
		}

		vector<double> getInitU(double h) {
			vector<double> U(N + 1);
			U[0] = 0;
			for (int i = 1; i <= N; ++i) {
				U[i] = phi(h*i);
			}
			return U;
		}

		vector<double> solveSystem(vector<double> Uprev, double tau, double h, int j) {
			vector<double> Unew(N + 1);
			Unew[0] = 0;
			for (int i = 1; i <= N; ++i) {
				double b = (tau*C(h*i, tau*j+tau)) / h;
				Unew[i] = (1 - b)*Uprev[i];
				if (i > 0) Unew[i] += b*Uprev[i - 1];
			}
			return Unew;
		}

	private:

		const int N, Nt;
		const fptr_d phi;
		const fptr_dd C, u;

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