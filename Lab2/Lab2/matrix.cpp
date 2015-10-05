#include "matrix.h"

#include <cassert>
#include <ostream>
#include <iostream>
#include <numeric>
#include <exception>

namespace Garbage {

	matrix::matrix(int m_rows, int m_columns)
		: m_rows(m_rows), m_columns(m_columns)
	{
		data.assign(m_rows, vector<_Type>(m_columns, 0));
	}

	matrix::matrix(const matrix& other)
		: m_rows(other.m_rows), m_columns(other.m_columns)
	{
		data.resize(m_rows);
		for (size_t i = 0; i < m_rows; ++i) {
			data[i].assign(other.data[i].begin(), other.data[i].end());
		}
	}

	matrix::~matrix() {
		for (size_t i = 0; i < data.size(); ++i) {
			data[i].clear();
		}
		data.clear();
	}

	inline
		size_t matrix::rows() const {
		return m_rows;
	}

	inline
		size_t matrix::columns() const {
		return m_columns;
	}

	vector<matrix::_Type> matrix::get_row(size_t i) const {
		if (i >= m_rows) {
			throw std::logic_error("\tget_row(size_t i) : invalid index.");
		}
		return data[i];
	}

	vector<matrix::_Type> matrix::get_column(size_t j) const {
		if (j >= m_columns) {
			throw std::logic_error("\tget_column(size_t i) : invalid index.");
		}
		vector<_Type> ret(m_rows);
		for (int i = 0; i < m_rows; ++i) {
			ret[i] = data[i][j];
		}
		return ret;
	}

	matrix matrix::transpos() const {
		matrix AT(m_columns, m_rows);
		for (int i = 0; i < AT.rows(); ++i) {
			for (int j = 0; j < AT.columns(); ++j) {
				AT[i][j] = data[j][i];
			}
		}
		return AT;
	}

	void matrix::operator=(const matrix& other) {
		m_rows = other.m_rows;
		m_columns = other.m_columns;
		data.resize(m_rows);
		for (size_t i = 0; i < m_rows; ++i) {
			data[i].assign(other.data[i].begin(), other.data[i].end());
		}
	}

	vector<matrix::_Type>& matrix::operator[](size_t i) {
		return data[i];
	}

	const vector<matrix::_Type>& matrix::operator[](size_t i) const {
		return data[i];
	}

	matrix matrix::operator+(const matrix& other) const {
		if (m_rows != other.m_rows || m_columns != other.m_columns) {
			throw std::logic_error("\tmatrix::operator+(const matrix& other) : matrices have a different size.");
		}
		matrix Res(m_rows, m_columns);
		for (size_t i = 0; i < m_rows; ++i) {
			for (size_t j = 0; j < m_columns; ++j) {
				Res[i][j] = data[i][j] + other.data[i][j];
			}
		}
		return Res;
	}

	matrix matrix::operator-(const matrix& other) const {
		if (m_rows != other.m_rows || m_columns != other.m_columns) {
			throw std::logic_error("\tmatrix::operator-(const matrix& other) : matrices have a different size.");
		}
		matrix Res(m_rows, m_columns);
		for (size_t i = 0; i < m_rows; ++i) {
			for (size_t j = 0; j < m_columns; ++j) {
				Res[i][j] = data[i][j] - other.data[i][j];
			}
		}
		return Res;
	}

	matrix matrix::operator*(const long double alpha) const {
		matrix Res(m_rows, m_columns);
		for (size_t i = 0; i < m_rows; ++i) {
			for (size_t j = 0; j < m_columns; ++j) {
				Res[i][j] = data[i][j] * alpha;
			}
		}
		return Res;
	}

	vector<matrix::_Type> matrix::operator*(vector<matrix::_Type>& v) const {
		int n = v.size();
		vector<_Type> res(n, 0);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				res[i] += data[i][j] * v[j];
			}
		}
		return vector<_Type>(res);
	}

	matrix matrix::operator*(const matrix& other) const {
		if (m_columns != other.m_rows) {
			throw std::logic_error("\tmatrix::operator*(const matrix& other) : matrices have an incompatible size.");
		}
		matrix Res(m_rows, other.m_columns);
		for (size_t i = 0; i < Res.m_rows; ++i) {
			for (size_t j = 0; j < Res.m_columns; ++j) {
				Res[i][j] = 0;
				for (size_t k = 0; k < m_columns; ++k) {
					Res[i][j] += data[i][k] * other.data[k][j];
				}
			}
		}
		return Res;
	}

	bool matrix::operator==(const matrix& other) {
		if (m_rows != other.m_rows || m_columns != other.m_columns) {
			return false;
		}
		const long double precission = 1e-5;
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_columns; ++j) {
				if (fabsl(data[i][j] - other.data[i][j]) > precission) {
					return false;
				}
			}
		}
		return true;
	}

	std::ostream& operator<<(std::ostream& os, const matrix& A) {
		for (int i = 0; i < A.m_rows; ++i) {
			for (int j = 0; j < A.m_columns; ++j) {
				os << (fabsl(A[i][j]) > 1e-5 ? A[i][j] : 0) << " ";
			}
			os << std::endl;
		}
		return os;
	}

	matrix matrix::LU_decomposition() const {
		if (m_rows != m_columns) {
			throw std::logic_error("\tmatrix::LU_decomposition() : required quadratic matrix.");
		}
		size_t n = m_rows;
		matrix L(n, n), U(n, n);
		for (int i = 0; i < n; ++i) {
			U[0][i] = data[0][i];
		}
		for (int i = 0; i < n; ++i) {
			L[i][0] = data[i][0] / U[0][0];
			L[i][i] = 1;
		}
		for (int i = 1; i < n; ++i) {
			for (int j = 1; j < n; ++j) {
				if (j < i) {
					L[i][j] = data[i][j];
					for (int k = 0; k < j; ++k) {
						L[i][j] -= L[i][k] * U[k][j];
					}
					L[i][j] /= U[j][j];
				}
			}
			for (int j = 1; j < n; ++j) {
				if (j >= i) {
					U[i][j] = data[i][j];
					for (int k = 0; k < i; ++k) {
						U[i][j] -= L[i][k] * U[k][j];
					}
				}
			}
		}
		for (int i = 0; i < n; ++i) {
			for (int j = i; j < n; ++j) {
				L[i][j] = U[i][j];
			}
		}
		return L;
	}

	void matrix::QR_decomposition(matrix& Q, matrix& R) const {
		if (m_rows != m_columns) {
			throw std::logic_error("\tmatrix::QR_decomposition() : required quadratic matrix.");
		}
		size_t n = m_rows;
		matrix E = matrix(n, n);
		Q = matrix(n, n);
		R = *this;
		for (int i = 0; i < n; ++i) {
			E[i][i] = 1;
			Q[i][i] = 1;
		}
		for (int i = 0; i < n - 1; ++i) {
			vector<_Type> cur_column = R.get_column(i);
			long double nonzero_tail = 0;
			for (int j = i + 1; j < n; ++j) {
				nonzero_tail += cur_column[j] * cur_column[j];
			}
			if (nonzero_tail) {
				long double cur_column_norm = 0;
				for (int j = i; j < n; ++j) {
					cur_column_norm += cur_column[j] * cur_column[j];
				}
				cur_column_norm = sqrtl(cur_column_norm);
				if (cur_column[i] >= 0) {
					cur_column[i] += cur_column_norm;
				}
				else {
					cur_column[i] -= cur_column_norm;
				}
				vector<matrix::_Type> u(n);
				for (int j = i; j < n; ++j) {
					u[j] = cur_column[j];
				}
				long double u_norm = 0;
				for (int j = i; j < n; ++j) {
					u_norm += u[j] * u[j];
				}
				u_norm = sqrtl(u_norm);
				for (int j = i; j < n; ++j) {
					u[j] /= u_norm;
				}
				vector<_Type> utR(n, 0);
				for (int j = 0; j < n; ++j) {
					for (int k = 0; k < n; ++k) {
						utR[j] += u[k] * R[k][j];
					}
				}
				for (int j = 0; j < n; ++j) {
					for (int k = 0; k < n; ++k) {
						R[j][k] -= 2 * u[j] * utR[k];
					}
				}
				vector<_Type> utQ(n, 0);
				for (int j = 0; j < n; ++j) {
					for (int k = 0; k < n; ++k) {
						utQ[j] += u[k] * Q[k][j];
					}
				}
				for (int j = 0; j < n; ++j) {
					for (int k = 0; k < n; ++k) {
						Q[j][k] -= 2 * u[j] * utQ[k];
					}
				}
			}
		}
		Q = Q.transpos();
	}

}