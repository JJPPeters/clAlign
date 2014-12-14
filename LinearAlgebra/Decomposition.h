#include "Matrix.h"
#include "Cholesky.h"
#include "LowerUpper.h"

template <typename T>
std::vector<T> lsSolver(const Matrix<T> &A, const std::vector<T> &b)
{
	assert(A.rows() == b.size());
	assert(A.rows() > A.cols());

	Cholesky<T> cho;
	cho.decompose(TransposeMultiply(A, A));
	if (cho.isSpd())
		return cho.solve(TransposeMultiply(A, b));
	else
		return luSolver(TransposeMultiply(A, A), TransposeMultiply(A, b));
}

template <typename T>
std::vector<T> luSolver(const Matrix<T> &A, const std::vector<T> &b)
{
	assert(A.rows() == b.size());

	LowerUpper<T> lu;
	lu.decompose(A);
	return lu.solve(b);
}


template <typename T>
Matrix<T> luSolver(const Matrix<T> &A, const Matrix<T> &B)
{
	assert(A.rows() == B.rows());

	LUD<Type> lu;
	lu.dec(A);
	return lu.solve(B);
}