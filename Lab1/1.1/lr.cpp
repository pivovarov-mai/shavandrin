#include <iostream>
#include <vector>

using namespace std;

int n;

void LU_decompose(vector<vector<float>> A, vector<vector<float>> &L, vector<vector<float>> &U)
{
	U = A;
	for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
			L[j][i] = U[j][i] / U[i][i];
	
	for (int k = 1; k < n; ++k)
	{
		for (int i = k - 1; i < n; ++i)
			for (int j = i; j < n; ++j)
				L[j][i] = U[j][i] / U[i][i];

		for (int i = k; i < n; ++i)
			for (int j = k-1; j < n; ++j)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}

}

void matrix_mult(vector<vector<float>> A, vector<vector<float>> B, vector<vector<float>> &R)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				R[i][j] += A[i][k] * B[k][j];
}

void print_matrix(vector<vector<float>> A)
{
	for (int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			cout << A[i][j] << " ";
		}
		cout << '\n';
	}
}

void print_vector(vector<float> A)
{
	cout << "(";
	for(int i = 0; i < n; ++i)
	{
		cout << A[i];
		if (i != n - 1)
			cout << ' ';
	}
	cout << ")\n";
}

float sumY(vector<vector<float>> L, vector<float> b, vector<float> y, int i)
{
	float res = 0;
	for (int j = 0; j < i; ++j)
		res += L[i][j] * y[j];
	return res;
}

float sumX(vector<vector<float>> U, vector<float> b, vector<float> x, int i)
{
	float res = 0;
	for (int j = n - 1; j > i; --j)
		res += U[i][j] * x[j];
	return res;
}

void solve(vector<vector<float>> L, vector<vector<float>> U, vector<float> b, vector<float> &y, vector<float> &x)
{
	for (int i = 0; i < n; ++i)
		y[i] = (b[i] -  sumY(L, b, y, i)) / L[i][i];

	for (int i = n - 1; i >= 0; --i)
		x[i] = (y[i] -  sumX(U, b, x, i)) / U[i][i];
}

int calc_det(vector<vector<float>> A)
{
	int res = 1;
	for (int i = 0; i < n; ++i)
		res *= A[i][i];
	return res;
}

void transpose(vector<vector<float>> &A)
{
	vector<vector<float>> AT(n, vector<float>(n, 0));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			AT[j][i] = A[i][j];
		}
	}
	A = AT;
}

int main()
{
	cout << "Enter n: ";
    cin >> n;
	vector<vector<float>> A(n, vector<float>(n, 0)), L(n, vector<float>(n, 0)), U(n, vector<float>(n, 0)), R(n, vector<float>(n, 0));
	vector<float> b(n, 0), y(n, 1), x(n, 1);
	cout << "Enter elements of matrix A:\n";
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			cin >> A[i][j];
	cout << "Enter elements of vector b:\n";
	for (int i = 0; i < n; ++i)
		cin >> b[i];
	LU_decompose(A,L,U);
	solve(L, U, b, y, x);
	cout << "x = ";
	print_vector(x);
	cout << "det A = ";
	int det = calc_det(U);
	cout << det << '\n';
	cout << "Inverse matrix: " << '\n';
	vector<vector<float>> rev_A(n, vector<float>(n, 0));
	for (int i = 0; i < n; ++i)
	{
		vector<float> b_k(n, 0), e_k(n, 0), z_k(n, 0);
		e_k[i] = 1;
		solve(L, U, e_k, z_k, b_k);
		rev_A[i] = b_k;
	}
	transpose(rev_A);
	print_matrix(rev_A);
	return 0;
}