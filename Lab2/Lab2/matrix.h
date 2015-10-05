#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <ostream>
#include <vector>
using std::vector;

namespace Garbage {

	class matrix {
	public:

		typedef long double _Type;

		matrix(int rows, int columns);
		matrix(const matrix& other);
		~matrix();

		size_t rows() const;
		size_t columns() const;

		vector<_Type> get_row(size_t i) const;
		vector<_Type> get_column(size_t j) const;

		matrix transpos() const;

		matrix LU_decomposition() const;
		void   QR_decomposition(matrix& Q, matrix& R) const;

		vector<_Type>& operator[](size_t i);
		const vector<_Type>& operator[](size_t i) const;

		void   operator=(const matrix& other);
		bool   operator==(const matrix& other);
		matrix operator+(const matrix& other) const;
		matrix operator-(const matrix& other) const;
		matrix operator*(const matrix& other) const;
		vector<_Type> operator*(vector<_Type>& v) const;
		matrix operator*(const long double alpha) const;

		friend std::ostream& operator<<(std::ostream& os, const matrix& A);

	private:

		size_t m_columns, m_rows;
		vector<vector<_Type>> data;
	};

}

#endif/*__MATRIX_H__*/