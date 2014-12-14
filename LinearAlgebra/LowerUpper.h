#include "Matrix.h"

template <typename T>
class LowerUpper
{
public:

	bool isNonSingular();

	std::vector<T> permuteCopy(const std::vector<T> &A, const std::vector<int> &piv);

	Matrix<T> permuteCopy(const Matrix<T> &A, const std::vector<int> &piv, int j0, int j1);

	void decompose(const Matrix<T> &A);

	std::vector<T> solve(const std::vector<T> &b);

	Matrix<T> solve(const Matrix<T> &B);

private:
	int m;
	int n;

	Matrix<T> LU;

	int pivsign;

	std::vector<int> piv;
};


template <typename T>
inline bool LowerUpper<T>::isNonSingular()
{
	for (int j = 0; j<n; ++j)
		if (abs(LU(j, j)) == 0)
			return false;

	return true;
}

template <typename T>
std::vector<T> LowerUpper<T>::permuteCopy(const std::vector<T> &A, const std::vector<int> &piv)
{
	int pivLength = piv.size();
	if (pivLength != A.size())
		return std::vector<T>();

	std::vector<T> x(pivLength);
	for (int i = 0; i<pivLength; ++i)
		x[i] = A[piv[i]];

	return x;
}

template <typename T>
Matrix<T> LowerUpper<T>::permuteCopy(const Matrix<T> &A, const std::vector<int> &piv, int j0, int j1)
{
	int pivLength = piv.size();
	Matrix<T> X(pivLength, j1 - j0 + 1);

	for (int i = 0; i<pivLength; ++i)
		for (int j = j0; j <= j1; ++j)
			X(i, j - j0) = A(piv[i], j);

	return X;
}

template <typename T>
void LowerUpper<T>::decompose(const Matrix<T> &A)
{
	m = A.rows();
	n = A.cols();
	piv = std::vector<int>(m);
	LU = A;

	// Use a "left-looking", dot-product, Crout/Doolittle algorithm.
	for (int i = 0; i<m; ++i)
		piv[i] = i;

	pivsign = 1;
	T *LUrowi = 0;
	std::vector<T> LUcolj(m);

	// outer loop
	for (int j = 0; j<n; ++j)
	{
		// Make a copy of the j-th column to localize references.
		for (int i = 0; i<m; ++i)
			LUcolj[i] = LU(i, j);

		// Apply previous transformations.
		for (int i = 0; i<m; ++i)
		{
			LUrowi = &LU(i, 0);

			// Most of the time is spent in the following dot product.
			int kmax = (i < j) ? i : j;
			T s = 0;

			for (int k = 0; k<kmax; ++k)
				s += LUrowi[k] * LUcolj[k];

			LUrowi[j] = LUcolj[i] -= s;
		}

		// Find pivot and exchange if necessary.
		int p = j;
		for (int i = j + 1; i<m; ++i)
			if (abs(LUcolj[i]) > abs(LUcolj[p]))
				p = i;

		if (p != j)
		{
			int k = 0;
			for (k = 0; k<n; ++k)
				std::swap(LU(p, k), LU(j, k));

			std::swap(piv[p], piv[j]);
			pivsign = -pivsign;
		}

		// compute multipliers
		if ((j < m) && (abs(LU(j, j)) != 0))
		for (int i = j + 1; i<m; ++i)
			LU(i, j) /= LU(j, j);
	}
}

template <typename T>
std::vector<T> LowerUpper<T>::solve(const std::vector<T> &b)
{
	// dimensions: A is mxn, X is nxk, B is mxk
	if (b.size() != m)
		return std::vector<T>();

	if (!isNonSingular())
		return std::vector<T>();

	std::vector<T> x = permuteCopy(b, piv);

	// solve L*Y = B(piv)
	for (int k = 0; k<n; ++k)
		for (int i = k + 1; i<n; ++i)
			x[i] -= x[k] * LU(i, k);

	// solve U*x = y;
	for (int k = n - 1; k >= 0; --k)
	{
		x[k] /= LU(k, k);
		for (int i = 0; i<k; ++i)
			x[i] -= x[k] * LU(i, k);
	}

	return x;
}

template <typename T>
Matrix<T> solve(const Matrix<T> &B)
{
	// dimensions: A is mxn, X is nxk, B is mxk
	if (B.rows() != m)
		return Matrix<T>();

	if (!isNonSingular())
		return Matrix<T>();

	// copy right hand side with pivoting
	int nx = B.cols();
	Matrix<T> X = permuteCopy(B, piv, 0, nx - 1);

	// solve L*Y = B(piv,:)
	for (int k = 0; k<n; ++k)
		for (int i = k + 1; i<n; ++i)
			for (int j = 0; j<nx; ++j)
				X(i, j) -= X(k, j) * LU(i, k);

	// solve U*X = Y;
	for (int k = n - 1; k >= 0; --k)
	{
		for (int j = 0; j<nx; ++j)
			X(k, j) /= LU(k, k);

		for (int i = 0; i<k; ++i)
			for (int j = 0; j<nx; ++j)
				X(i, j) -= X(k, j) * LU(i, k);
	}

	return X;
}