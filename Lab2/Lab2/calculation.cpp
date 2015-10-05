#include "calculation.h"

#include <iostream>

namespace Garbage {

	vector<matrix::_Type> calculation::SLAE::back_substitution(const matrix& U, const vector<matrix::_Type>& b) {
		vector<matrix::_Type> x(U.columns());
		for (int i = x.size() - 1; i >= 0; --i) {
			x[i] = b[i];
			for (int j = i + 1; j < U.columns(); ++j) {
				x[i] -= U[i][j] * x[j];
			}
			x[i] /= U[i][i];
		}
		return x;
	}

	vector<matrix::_Type> calculation::SLAE::straight_substitution(const matrix& L, const vector<matrix::_Type>& b) {
		vector<matrix::_Type> x(L.columns());
		for (int i = 0; i < x.size(); ++i) {
			x[i] = b[i];
			for (int j = 0; j < i; ++j) {
				x[i] -= L[i][j] * x[j];
			}
			x[i] /= L[i][i];
		}
		return x;
	}

	matrix calculation::SLAE::LU_method(const matrix& A, const matrix& b) {
		matrix LU = A.LU_decomposition();
		matrix L(LU.rows(), LU.columns()), U(LU.rows(), LU.columns());
		for (int i = 0; i < LU.rows(); ++i) {
			for (int j = 0; j < i; ++j) {
				L[i][j] = LU[i][j];
			}
			L[i][i] = 1;
		}
		for (int i = 0; i < LU.columns(); ++i) {
			for (int j = i; j < LU.columns(); ++j) {
				U[i][j] = LU[i][j];
			}
		}
		//LU.~~matrix();
		matrix X(A.columns(), b.columns());
		for (int i = 0; i < b.columns(); ++i) {
			vector<matrix::_Type> cur_column = b.get_column(i);
			vector<matrix::_Type> y = straight_substitution(L, cur_column);
			vector<matrix::_Type> x = back_substitution(U, y);
			for (int j = 0; j < X.rows(); ++j) {
				X[j][i] = x[j];
			}
		}
		return X;
	}

	matrix calculation::SLAE::QR_method(const matrix& A, const matrix& b) {
		//if (A.rows != A.columns) {
		//
		//}
		int n = A.rows();
		matrix Q(n, n), R(n, n);
		A.QR_decomposition(Q, R);
		matrix X(A.columns(), b.columns());
		matrix QTb = Q.transpos() * b;
		for (int i = 0; i < b.columns(); ++i) {
			/*vector<matrix::_Type> cur_column = b.get_column(i);
			vector<matrix::_Type> y(n, 0);
			for (int j = 0; j < n; ++j) {
			for (int k = 0; k < n; ++k) {
			y[j] -= Q[k][j] * cur_column[k];
			}
			}
			vector<matrix::_Type> x = back_substitution(R, y);*/
			vector<matrix::_Type> cur_column = QTb.get_column(i);
			vector<matrix::_Type> x = back_substitution(R, cur_column);
			for (int j = 0; j < X.rows(); ++j) {
				X[j][i] = x[j];
			}
		}
		return X;
	}

	matrix calculation::inverse_matrix::LU_method(const matrix& A) {
		//if (A.rows != A.columns) {
		//
		//}
		matrix E(A.rows(), A.columns());
		for (int i = 0; i < A.rows(); ++i) {
			E[i][i] = 1;
		}
		matrix A_inverse = SLAE::LU_method(A, E);
		return A_inverse;
	}

	matrix calculation::inverse_matrix::QR_method(const matrix& A) {
		//if (A.rows != A.columns) {
		//
		//}
		matrix E(A.rows(), A.columns());
		for (int i = 0; i < A.rows(); ++i) {
			E[i][i] = 1;
		}
		matrix A_inverse = SLAE::QR_method(A, E);
		return A_inverse;
	}

}