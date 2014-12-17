#include "Matrix.h"
#include<complex>

const float EPS = 2.220446049250313e-016;

template <typename T>
class SVD
{
public:
	void decompose(const Matrix<T> &A);

	int rank();

private:
	Matrix<T> U;
	Matrix<T> V;
	std::vector<T> S;

	void doDecomposition(Matrix<T> &B, Matrix<T> &U, std::vector<T> &S, Matrix<T> &V);
};


template <typename T>
class SVD<std::complex<T>>
{
public:
	void decompose(const Matrix<std::complex<T>> &A);

	int rank();

private:
	Matrix<std::complex<T>> U;
	Matrix<std::complex<T>> V;
	std::vector<T> S;

	void doDecomposition(Matrix<std::complex<T>> &B, Matrix<std::complex<T>> &U, std::vector<T> &S, Matrix<std::complex<T>> &V);
};

template <typename T>
void SVD<T>::decompose(const Matrix<T> &A)
{
	int m = A.rows(),
		n = A.cols(),
		p = std::min(m, n);

	U = Matrix<T>(m, p);
	V = Matrix<T>(n, p);
	S = std::vector<T>(p);
	if (m >= n)
	{
		Matrix<T> B(A);
		doDecomposition(B, U, S, V);
	}
	else
	{
		Matrix<T> B(Transpose(A));
		doDecomposition(B, V, S, U);
	}
}

template <typename T>
void SVD<std::complex<T>>::decompose(const Matrix<std::complex<T>> &A)
{
	int m = A.rows(),
		n = A.cols(),
		p = std::min(m, n);

	U = Matrix<std::complex<T>>(m, p);
	V = Matrix<std::complex<T>>(n, p);
	S = std::vector<T>(p);
	if (m >= n)
	{
		Matrix<std::complex<T>> B(A);
		doDecomposition(B, U, S, V);
	}
	else
	{
		Matrix<std::complex<T>> B(Transpose(A));
		doDecomposition(B, V, S, U);
	}
}

template <typename T>
int SVD<T>::rank()
{
	int N = S.size();
	float tol = N * S[0] * EPS;
	int r = 0;

	for (int i = 0; i<N; ++i)
		if (S[i] > tol)
			r++;

	return r;
}

template <typename T>
int SVD<std::complex<T>>::rank()
{
	int N = S.size();
	float tol = N * S[0] * EPS;
	int r = 0;

	for (int i = 0; i<N; ++i)
		if (S[i] > tol)
			r++;

	return r;
}

template <typename T>
void SVD<T>::doDecomposition(Matrix<T> &B, Matrix<T> &U, std::vector<T> &S, Matrix<T> &V)
{
	int m = B.rows(),
		n = B.cols();

	std::vector<T> e(n);
	std::vector<T> work(m);

	// boolean
	int wantu = 1;
	int wantv = 1;

	// Reduce A to bidiagonal form, storing the diagonal elements
	// in s and the super-diagonal elements in e.
	int nct = std::min(m - 1, n);
	int nrt = std::max(0, n - 2);
	int i = 0,
		j = 0,
		k = 0;

	for (k = 0; k<std::max(nct, nrt); ++k)
	{
		if (k < nct)
		{
			// Compute the transformation for the k-th column and
			// place the k-th diagonal in s[k].
			// Compute 2-norm of k-th column without under/overflow.
			S[k] = 0;
			for (i = k; i<m; ++i)
				S[k] = hypot(S[k], B(i, k));

			if (S[k] != 0)
			{
				if (B(k, k) < 0)
					S[k] = -S[k];

				for (i = k; i<m; ++i)
					B(i, k) /= S[k];
				B(k, k) += 1;
			}
			S[k] = -S[k];
		}

		for (j = k + 1; j<n; ++j)
		{
			if ((k < nct) && (S[k] != 0))
			{
				// apply the transformation
				T t = 0;
				for (i = k; i<m; ++i)
					t += B(i, k) * B(i, j);

				t = -t / B(k, k);
				for (i = k; i<m; ++i)
					B(i, j) += t*B(i, k);
			}
			e[j] = B(k, j);
		}

		// Place the transformation in U for subsequent back
		// multiplication.
		if (wantu & (k < nct))
			for (i = k; i<m; ++i)
				U(i, k) = B(i, k);

		if (k < nrt)
		{
			// Compute the k-th row transformation and place the
			// k-th super-diagonal in e[k].
			// Compute 2-norm without under/overflow.
			e[k] = 0;
			for (i = k + 1; i<n; ++i)
				e[k] = hypot(e[k], e[i]);

			if (e[k] != 0)
			{
				if (e[k + 1] < 0)
					e[k] = -e[k];

				for (i = k + 1; i<n; ++i)
					e[i] /= e[k];
				e[k + 1] += 1;
			}
			e[k] = -e[k];

			if ((k + 1 < m) && (e[k] != 0))
			{
				// apply the transformation
				for (i = k + 1; i<m; ++i)
					work[i] = 0;

				for (j = k + 1; j<n; ++j)
					for (i = k + 1; i<m; ++i)
						work[i] += e[j] * B(i, j);

				for (j = k + 1; j<n; ++j)
				{
					T t = -e[j] / e[k + 1];
					for (i = k + 1; i<m; ++i)
						B(i, j) += t * work[i];
				}
			}

			// Place the transformation in V for subsequent
			// back multiplication.
			if (wantv)
				for (i = k + 1; i<n; ++i)
					V(i, k) = e[i];
		}
	}

	// Set up the final bidiagonal matrix or order p.
	int p = n;

	if (nct < n)
		S[nct] = B(nct, nct);
	if (m < p)
		S[p - 1] = 0;

	if (nrt + 1 < p)
		e[nrt] = B(nrt, p - 1);
	e[p - 1] = 0;

	// if required, generate U
	if (wantu)
	{
		for (j = nct; j<n; ++j)
		{
			for (i = 0; i<m; ++i)
				U(i, j) = 0;
			U(j, j) = 1;
		}

		for (k = nct - 1; k >= 0; --k)
			if (S[k] != 0)
			{
				for (j = k + 1; j<n; ++j)
				{
					T t = 0;
					for (i = k; i<m; ++i)
						t += U(i, k) * U(i, j);
					t = -t / U(k, k);

					for (i = k; i<m; ++i)
						U(i, j) += t * U(i, k);
				}

				for (i = k; i<m; ++i)
					U(i, k) = -U(i, k);
				U(k, k) = 1 + U(k, k);

				for (i = 0; i<k - 1; ++i)
					U(i, k) = 0;
			}
			else
			{
				for (i = 0; i<m; ++i)
					U(i, k) = 0;
				U(k, k) = 1;
			}
	}

	// if required, generate V
	if (wantv)
		for (k = n - 1; k >= 0; --k)
		{
			if ((k < nrt) && (e[k] != 0))
				for (j = k + 1; j<n; ++j)
				{
					T t = 0;
					for (i = k + 1; i<n; ++i)
						t += V(i, k) * V(i, j);
					t = -t / V(k + 1, k);

					for (i = k + 1; i<n; ++i)
						V(i, j) += t * V(i, k);
				}

			for (i = 0; i<n; ++i)
				V(i, k) = 0;
			V(k, k) = 1;
		}

	// main iteration loop for the singular values
	int pp = p - 1;
	int iter = 0;
	float eps = pow(2.0, -52.0);

	while (p > 0)
	{
		int k = 0;
		int kase = 0;

		// Here is where a test for too many iterations would go.
		// This section of the program inspects for negligible
		// elements in the s and e arrays. On completion the
		// variables kase and k are set as follows.
		// kase = 1     if s(p) and e[k-1] are negligible and k<p
		// kase = 2     if s(k) is negligible and k<p
		// kase = 3     if e[k-1] is negligible, k<p, and
		//				s(k), ..., s(p) are not negligible
		// kase = 4     if e(p-1) is negligible (convergence).
		for (k = p - 2; k >= -1; --k)
		{
			if (k == -1)
				break;

			if (abs(e[k]) <= eps*(abs(S[k]) + abs(S[k + 1])))
			{
				e[k] = 0;
				break;
			}
		}

		if (k == p - 2)
			kase = 4;
		else
		{
			int ks;
			for (ks = p - 1; ks >= k; --ks)
			{
				if (ks == k)
					break;

				T t = ((ks != p) ? abs(e[ks]) : 0) +
					((ks != k + 1) ? abs(e[ks - 1]) : 0);

				if (abs(S[ks]) <= eps*t)
				{
					S[ks] = 0;
					break;
				}
			}

			if (ks == k)
				kase = 3;
			else if (ks == p - 1)
				kase = 1;
			else
			{
				kase = 2;
				k = ks;
			}
		}
		k++;

		// Perform the task indicated by kase.
		switch (kase)
		{
			// deflate negligible s(p)
		case 1:
		{
			T f = e[p - 2];
			e[p - 2] = 0;

			for (j = p - 2; j >= k; --j)
			{
				T t = hypot(S[j], f);
				T cs = S[j] / t;
				T sn = f / t;
				S[j] = t;

				if (j != k)
				{
					f = -sn * e[j - 1];
					e[j - 1] = cs * e[j - 1];
				}

				if (wantv)
					for (i = 0; i<n; ++i)
					{
						t = cs*V(i, j) + sn*V(i, p - 1);
						V(i, p - 1) = -sn*V(i, j) + cs*V(i, p - 1);
						V(i, j) = t;
					}
			}
		}
		break;

		// split at negligible s(k)
		case 2:
		{
			T f = e[k - 1];
			e[k - 1] = 0;

			for (j = k; j<p; ++j)
			{
				T t = hypot(S[j], f);
				T cs = S[j] / t;
				T sn = f / t;
				S[j] = t;
				f = -sn * e[j];
				e[j] = cs * e[j];

				if (wantu)
					for (i = 0; i<m; ++i)
					{
						t = cs*U(i, j) + sn*U(i, k - 1);
						U(i, k - 1) = -sn*U(i, j) + cs*U(i, k - 1);
						U(i, j) = t;
					}
			}
		}
		break;

		// perform one qr step
		case 3:
		{
			// calculate the shift
			T scale = std::max(std::max(std::max(std::max(
				abs(S[p - 1]), abs(S[p - 2])), abs(e[p - 2])),
				abs(S[k])), abs(e[k]));
			T sp = S[p - 1] / scale;
			T spm1 = S[p - 2] / scale;
			T epm1 = e[p - 2] / scale;
			T sk = S[k] / scale;
			T ek = e[k] / scale;
			T b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1) / 2.0;
			T c = (sp*epm1) * (sp*epm1);
			T shift = 0;

			if ((b != 0) || (c != 0))
			{
				shift = sqrt(b*b + c);
				if (b < 0)
					shift = -shift;
				shift = c / (b + shift);
			}
			T f = (sk + sp)*(sk - sp) + shift;
			T g = sk * ek;

			// chase zeros
			for (j = k; j<p - 1; ++j)
			{
				T t = hypot(f, g);
				T cs = f / t;
				T sn = g / t;
				if (j != k)
					e[j - 1] = t;

				f = cs*S[j] + sn*e[j];
				e[j] = cs*e[j] - sn*S[j];
				g = sn * S[j + 1];
				S[j + 1] = cs * S[j + 1];

				if (wantv)
					for (i = 0; i<n; ++i)
					{
						t = cs*V(i, j) + sn*V(i, j + 1);
						V(i, j + 1) = -sn*V(i, j) + cs*V(i, j + 1);
						V(i, j) = t;
					}

				t = hypot(f, g);
				cs = f / t;
				sn = g / t;
				S[j] = t;
				f = cs*e[j] + sn*S[j + 1];
				S[j + 1] = -sn*e[j] + cs*S[j + 1];
				g = sn * e[j + 1];
				e[j + 1] = cs * e[j + 1];

				if (wantu && (j < m - 1))
					for (i = 0; i<m; ++i)
					{
						t = cs*U(i, j) + sn*U(i, j + 1);
						U(i, j + 1) = -sn*U(i, j) + cs*U(i, j + 1);
						U(i, j) = t;
					}
			}
			e[p - 2] = f;
			iter = iter + 1;
		}
		break;

		// convergence
		case 4:
		{
			// Make the singular values positive.
			if (S[k] <= 0)
			{
				S[k] = (S[k] < 0) ? -S[k] : 0;
				if (wantv)
					for (i = 0; i <= pp; ++i)
						V(i, k) = -V(i, k);
			}

			// Order the singular values.
			while (k < pp)
			{
				if (S[k] >= S[k + 1])
					break;

				T t = S[k];
				S[k] = S[k + 1];
				S[k + 1] = t;

				if (wantv && (k < n - 1))
					for (i = 0; i<n; ++i)
						std::swap(V(i, k), V(i, k + 1));

				if (wantu && (k < m - 1))
					for (i = 0; i<m; ++i)
						std::swap(U(i, k), U(i, k + 1));
				k++;
			}
			iter = 0;
			p--;
		}
		break;
		}
	}
}


template <typename T>
void SVD<std::complex<T>>::doDecomposition(Matrix<std::complex<T>> &B, Matrix<std::complex<T>> &U, std::vector<T> &S, Matrix<std::complex<T>> &V)
{
	int k, k1,
		L, L1,
		nM1,
		m = B.rows(),
		n = B.cols();

	T    tol, eta;
	T    w, x, y, z;
	T    cs, sn, f, g, h;
	std::complex<T>   q;

	std::vector<T>    b(m), c(m), t(m);

	L = 0;
	nM1 = n - 1;
	eta = 2.8E-16;      // the relative machine precision

	// Householder Reduction
	c[0] = 0;
	k = 0;
	while (1)
	{
		k1 = k + 1;

		// elimination of a(i, k), i = k, ..., m-1
		z = 0;
		for (int i = k; i<m; ++i)
			z += std::norm(B(i, k));
		b[k] = 0;

		if (z > EPS)
		{
			z = sqrt(z);
			b[k] = z;
			w = abs(B(k, k));
			q = 1;
			if (w != T(0))
				q = B(k, k) / w;

			B(k, k) = q * (z + w);

			if (k != n - 1)
				for (int j = k1; j<n; ++j)
				{
					q = 0;
					for (int i = k; i<m; ++i)
						q += conj(B(i, k)) * B(i, j);
					q /= z * (z + w);

					for (int i = k; i<m; ++i)
						B(i, j) -= q * B(i, k);
				}

			// phase transformation
			q = -conj(B(k, k)) / abs(B(k, k));
			for (int j = k1; j<n; ++j)
				B(k, j) *= q;
		}

		// elimination of a(k, j), j = k+2, ..., n-1
		if (k == nM1)
			break;
		z = 0;
		for (int j = k1; j<n; ++j)
			z += norm(B(k, j));
		c[k1] = 0;

		if (z > EPS)
		{
			z = sqrt(z);
			c[k1] = z;
			w = abs(B(k, k1));
			q = 1;
			if (w != T(0))
				q = B(k, k1) / w;

			B(k, k1) = q * (z + w);
			for (int i = k1; i<m; ++i)
			{
				q = 0;
				for (int j = k1; j<n; ++j)
					q += conj(B(k, j)) * B(i, j);

				q /= z * (z + w);
				for (int j = k1; j<n; ++j)
					B(i, j) -= q * B(k, j);
			}

			// phase transformation
			q = -conj(B(k, k1)) / abs(B(k, k1));
			for (int i = k1; i<m; ++i)
				B(i, k1) *= q;
		}
		k = k1;
	}

	// tolerance for negligible elements
	tol = 0;
	for (k = 0; k<n; ++k)
	{
		S[k] = b[k];
		t[k] = c[k];
		if (S[k] + t[k] > tol)
			tol = S[k] + t[k];
	}
	tol *= eta;

	// initialization fo U and V
	for (int j = 0; j<n; ++j)
	{
		for (int i = 0; i<m; ++i)
			U(i, j) = 0;
		U(j, j) = 1;

		for (int i = 0; i<n; ++i)
			V(i, j) = 0;
		V(j, j) = 1;
	}

	// QR diagonallization
	for (k = nM1; k >= 0; --k)
	{
		// test for split
		while (1)
		{
			for (L = k; L >= 0; --L)
			{
				if (abs(t[L]) <= tol)
					goto Test;
				if (abs(S[L - 1]) <= tol)
					break;
			}

			// cancellation of E(L)
			cs = 0;
			sn = 1;
			L1 = L - 1;
			for (int i = L; i <= k; ++i)
			{
				f = sn * t[i];
				t[i] *= cs;
				if (abs(f) <= tol)
					goto Test;

				h = S[i];
				w = sqrt(f*f + h*h);
				S[i] = w;
				cs = h / w;
				sn = -f / w;

				for (int j = 0; j<n; ++j)
				{
					x = real(U(j, L1));
					y = real(U(j, i));
					U(j, L1) = x*cs + y*sn;
					U(j, i) = y*cs - x*sn;
				}
			}

			// test for convergence
		Test:	w = S[k];
			if (L == k)
				break;

			// origin shift
			x = S[L];
			y = S[k - 1];
			g = t[k - 1];
			h = t[k];
			f = ((y - w)*(y + w) + (g - h)*(g + h)) / (2 * h*y);
			g = sqrt(f*f + 1);
			if (f < T(0))
				g = -g;
			f = ((x - w)*(x + w) + h*(y / (f + g) - h)) / x;

			// QR step
			cs = 1;
			sn = 1;
			L1 = L + 1;
			for (int i = L1; i <= k; ++i)
			{
				g = t[i];
				y = S[i];
				h = sn * g;
				g = cs * g;
				w = sqrt(h*h + f*f);
				t[i - 1] = w;
				cs = f / w;
				sn = h / w;
				f = x * cs + g * sn;
				g = g * cs - x * sn;
				h = y * sn;
				y = y * cs;

				for (int j = 0; j<n; ++j)
				{
					x = real(V(j, i - 1));
					w = real(V(j, i));
					V(j, i - 1) = x*cs + w*sn;
					V(j, i) = w*cs - x*sn;
				}

				w = sqrt(h*h + f*f);
				S[i - 1] = w;
				cs = f / w;
				sn = h / w;
				f = cs * g + sn * y;
				x = cs * y - sn * g;

				for (int j = 0; j<n; ++j)
				{
					y = real(U(j, i - 1));
					w = real(U(j, i));
					U(j, i - 1) = y*cs + w*sn;
					U(j, i) = w*cs - y*sn;
				}
			}
			t[L] = 0;
			t[k] = f;
			S[k] = x;
		}

		// convergence
		if (w >= T(0))
			continue;
		S[k] = -w;

		for (int j = 0; j<n; ++j)
			V(j, k) = -V(j, k);
	}

	// sort dingular values
	for (k = 0; k<n; ++k)
	{
		g = -1.0;
		int j = k;
		for (int i = k; i<n; ++i)
		{
			if (S[i] <= g)
				continue;
			g = S[i];
			j = i;
		}

		if (j == k)
			continue;
		S[j] = S[k];
		S[k] = g;

		for (int i = 0; i<n; ++i)
		{
			q = V(i, j);
			V(i, j) = V(i, k);
			V(i, k) = q;
		}

		for (int i = 0; i<n; ++i)
		{
			q = U(i, j);
			U(i, j) = U(i, k);
			U(i, k) = q;
		}
	}

	// back transformation
	for (k = nM1; k >= 0; --k)
	{
		if (b[k] == T(0))
			continue;
		q = -B(k, k) / abs(B(k, k));

		for (int j = 0; j<n; ++j)
			U(k, j) *= q;

		for (int j = 0; j<n; ++j)
		{
			q = 0;
			for (int i = k; i<m; ++i)
				q += conj(B(i, k)) * U(i, j);
			q /= abs(B(k, k)) * b[k];

			for (int i = k; i<m; ++i)
				U(i, j) -= q * B(i, k);
		}
	}

	if (n > 1)
	{
		for (k = n - 2; k >= 0; --k)
		{
			k1 = k + 1;
			if (c[k1] == T(0))
				continue;
			q = -conj(B(k, k1)) / abs(B(k, k1));

			for (int j = 0; j<n; ++j)
				V(k1, j) *= q;

			for (int j = 0; j<n; ++j)
			{
				q = 0;
				for (int i = k1; i<n; ++i)
					q += B(k, i) * V(i, j);
				q /= abs(B(k, k1)) * c[k1];

				for (int i = k1; i<n; ++i)
					V(i, j) -= q * conj(B(k, i));
			}
		}
	}
}