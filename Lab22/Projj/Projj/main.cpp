#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

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

double phi(double x) {
	return sin(2 * M_PI*x) + 2.5*M_PI*x;
}

double C(double x, double t) {
	return (sin(2 * M_PI*x)*sin(M_PI*t)) / (2 * cos(2 * M_PI*x)*cos(M_PI*t) + 2.5*M_PI);
}

double u(double x, double t) {
	return sin(2 * M_PI*x)*cos(M_PI*t) + 2.5*M_PI*x;
}

/* b = (tau * C^j_k) / h
матрица n-1 x n-1
1-b 000000000000000
b 1-b 0000000000000
0 b 1-b 00000000000
00 b 1-b 0000000000
u = (u1, ... ,uN) */

int main() {
	freopen("output.txt", "w", stdout);
	int N = 10;
	int Nt = 100;
	Calculation calc(phi, C, u, N, Nt);
	calc.calculate();
	return 0;
}