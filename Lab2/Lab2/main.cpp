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