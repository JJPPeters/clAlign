#include "Matrix.h"
#include<complex>

template <typename T>
class Cholesky
{
public:

	bool isSpd() { return spd; }

	void decompose(const Matrix<T> &A);

	std::vector<T> solve(const std::vector<T> &b);

	Matrix<T> solve(const Matrix<T> &B);

private:
	bool spd;
	Matrix<T> L;
};

// Specialised Template for complex types
template <class T>
class Cholesky < std::complex<T> >
{
public:
	bool isSpd() { return spd; }

	void decompose(const Matrix<std::complex<T>> &A);

	std::vector<std::complex<T>> solve(const std::vector<std::complex<T>> &b);

	Matrix<std::complex<T>> solve(const Matrix<std::complex<T>> &B);

private:
	bool spd;
	Matrix<std::complex<T>> L;
};

template <typename T>
void Cholesky<T>::decompose(const Matrix<T> &A)
{
	int m = A.rows();
	int n = A.cols();

	spd = (m == n);
	if (!spd)
		return;

	L = Matrix<T>(n, n);

	for (int j = 0; j<m; ++j) // should be j < m ? originally < n
	{
		T d = 0;
		for (int k = 0; k<j; ++k)
		{
			T s = 0;
			for (int i = 0; i<k; ++i)
				s += L(k, i) * L(j, i);

			L(j, k) = s = (A(j, k) - s) / L(k, k);
			d = d + s*s;
			spd = spd && (A(k, j) == A(j, k));
		}

		d = A(j, j) - d;
		spd = spd && (d > 0);

		L(j, j) = sqrt(d > 0 ? d : 0);
		for (int k = j + 1; k<m; ++k) // should be k < m ? originally < n
			L(j, k) = 0;
	}
}

template <typename T>
void Cholesky<std::complex<T>>::decompose(const Matrix<std::complex<T>> &A)
{
	int m = A.rows();
	int n = A.cols();

	spd = (m == n);
	if (!spd)
		return;

	L = Matrix<std::complex<T>>(n, n);

	for (int j = 0; j<m; ++j)
	{
		spd = spd && (std::imag(A(j, j)) == 0);
		T d = 0;
		for (int k = 0; k<j; ++k)
		{
			std::complex<T> s = 0;
			for (int i = 0; i<k; ++i)
				s += L(k, i) * std::conj(L(j, i));

			L(j, k) = s = (A(j, k) - s) / L(k, k);
			d = d + std::norm(s);
			spd = spd && (A(k, j) == std::conj(A(j, k)));
		}

		d = std::real(A(j, j)) - d;
		spd = spd && (d > 0);

		L(j, j) = sqrt(d > 0 ? d : 0);
		for (int k = j + 1; k<m; ++k)
			L(j, k) = 0;
	}
}

template <typename T>
std::vector<T> Cholesky<T>::solve(const std::vector<T> &b)
{
	int n = L.rows();
	if (b.size() != n)
		return std::vector<T>();

	std::vector<T> x = b;

	// solve L*y = b
	for (int k = 0; k<n; ++k)
	{
		for (int i = 0; i<k; ++i)
			x[k] -= x[i] * L(k, i);

		x[k] /= L(k, k);
	}

	// solve L'*x = y
	for (int k = n - 1; k >= 0; --k)
	{
		for (int i = k + 1; i<n; ++i)
			x[k] -= x[i] * L(i, k);

		x[k] /= L(k, k);
	}

	return x;
}

template <typename T>
std::vector<std::complex<T>> Cholesky<std::complex<T>>::solve(const std::vector<std::complex<T>> &b)
{
	int n = L.rows();
	if (b.size() != n)
		return std::vector<std::complex<T>>();

	std::vector<std::complex<T>> x = b;

	// solve L*y = b
	for (int k = 0; k<n; ++k)
	{
		for (int i = 0; i<k; ++i)
			x[k] -= x[i] * L(k, i);

		x[k] /= L(k, k);
	}

	// solve L'*x = y
	for (int k = n - 1; k >= 0; --k)
	{
		for (int i = k + 1; i<n; ++i)
			x[k] -= x[i] * std::conj(L(i, k));

		x[k] /= L(k, k);
	}

	return x;
}

template <typename T>
Matrix<T> Cholesky<T>::solve(const Matrix<T> &B)
{
	int n = L.rows();
	if (B.rows() != n)
		return Matrix<T>();

	Matrix<T> X = B;
	int nx = B.cols();

	// solve L*Y = B
	for (int j = 0; j<nx; ++j)
		for (int k = 0; k<n; ++k)
		{
			for (int i = 0; i<k; ++i)
				X(k, j) -= X(i, j) * L(k, i);

			X(k, j) /= L(k, k);
		}

	// solve L'*X = Y
	for (int j = 0; j<nx; ++j)
		for (int k = n - 1; k >= 0; --k)
		{
			for (int i = k + 1; i<n; ++i)
				X(k, j) -= X(i, j) * L(i, k);

			X(k, j) /= L(k, k);
		}

	return X;
}

template <typename T>
Matrix<std::complex<T>> Cholesky<std::complex<T>>::solve(const Matrix<std::complex<T>> &B)
{
	int n = L.rows();
	if (B.rows() != n)
		return Matrix<std::complex<T>>();

	Matrix<std::complex<T>> X = B;
	int nx = B.cols();

	// solve L*Y = B
	for (int j = 0; j<nx; ++j)
		for (int k = 0; k<n; ++k)
		{
			for (int i = 0; i<k; ++i)
				X(k, j) -= X(i, j) * L(k, i);

			X(k, j) /= L(k, k);
		}

	// solve L'*X = Y
	for (int j = 0; j<nx; ++j)
		for (int k = n - 1; k >= 0; --k)
		{
			for (int i = k + 1; i<n; ++i)
				X(k, j) -= X(i, j) * std::conj(L(i, k));

			X(k, j) /= L(k, k);
		}

	return X;
}