#ifndef __CALCULATION_H__
#define __CALCULATION_H__

#include "matrix.h"

namespace Garbage {

	class calculation {
	public:

		class SLAE {
		public:

			static vector<matrix::_Type> back_substitution(const matrix& U, const vector<matrix::_Type>& b);
			static vector<matrix::_Type> straight_substitution(const matrix& L, const vector<matrix::_Type>& b);
			static matrix LU_method(const matrix& A, const matrix& b);
			static matrix QR_method(const matrix& A, const matrix& b);
		};

		class inverse_matrix {
		public:

			static matrix LU_method(const matrix& A);
			static matrix QR_method(const matrix& A);
		};

	};

}

#endif/*__CALCULATION_H__*/