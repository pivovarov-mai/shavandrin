#include <iostream>
#include <vector>

using namespace std;

const int SIZE_LIMIT = 1e3+13; // Ограничение на размер матрицы
int N; //Размер матрицы

vector <vector <double>> A (SIZE_LIMIT, vector <double>(SIZE_LIMIT, 0.0));
vector <vector <double>> A1 (SIZE_LIMIT, vector <double>(SIZE_LIMIT, 0.0));
vector <vector <double>> Y1 (SIZE_LIMIT, vector <double>(SIZE_LIMIT, 0.0));
vector <vector <double>> L (SIZE_LIMIT, vector <double>(SIZE_LIMIT, 0.0));
vector <vector <double>> U (SIZE_LIMIT, vector <double>(SIZE_LIMIT, 0.0));
vector <vector <double>> LU (SIZE_LIMIT, vector <double>(SIZE_LIMIT, 0.0));
vector <double> b (SIZE_LIMIT);
vector <double> x (SIZE_LIMIT);
vector <double> y (SIZE_LIMIT);

void ReadData()
{
    cout<<"? matrix size\n";
    cin>>N;

    cout<<"? input matrix\n";
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            cin>>A[i][j];

    cout<<"? answer matrix\n";
    for(int i = 0; i < N; i++)
        cin>>b[i];
}

void LU_d()
{
    for (int k = 0; k < N-1; k++) {
        for (int i = k; i < N; i++)
            for (int j = i; j < N; j++)
                L[j][i] = U[j][i] / U[i][i];

        for (int i = k + 1; i < N; i++)
            for (int j = k; j < N; j++)
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
    }
}

void Check()
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                LU[i][j] += L[i][k] * U[k][j];
}

void Solve()
{
    for (int i = 0; i < N; i++) {
        double S = 0;
        for (int j = 0; j < i; j++)
            S += L[i][j] * y[j];
        y[i] = b[i] - S;
    }

    for (int i = N - 1; i >= 0; i--) {
        double S = 0;
        for (int j = N - 1; j > i; j--)
            S += U[i][j] * x[j];
        x[i] = (y[i] - S) / U[i][i];
    }
}

void Inverse()
{
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < N; i++) {
            double S = 0;
            for (int j = 0; j < i; j++)
                S += L[i][j] * Y1[j][k];
            Y1[i][k] = ((i == k) ?  1 - S : -S);
        }
    }

    for (int k = 0; k < N; k++) {
        for (int i = N - 1; i >= 0; i--) {
            double S = 0;
            for (int j = N - 1; j > i; j--)
                S += U[i][j] * A1[j][k];
            A1[i][k] = (Y1[i][k] - S) / U[i][i];
        }
    }
}

double Det(vector <vector<double>> v) {
    double p = 1.0;
    for (int i = 0; i < N; i++)
        p *= v[i][i];
    return p;
}

void Print_matrix(vector <vector<double>> &v) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << v[i][j] << "\t\t";
        cout << "\n";
    }
}

void WriteData()
{
    cout << "\n\t\t\tMatrix L\n";
	cout << "-------------------------------------------------------------------\n";
	Print_matrix(L);
	cout << "\n\n\t\t\tMatrix U\n";
	cout << "-------------------------------------------------------------------\n";
	Print_matrix(U);
	cout << "\n\n\t\t\tMatrix L*U\n";
	cout << "-------------------------------------------------------------------\n";
	Print_matrix(LU);
	cout << "\n\n\t\t\tMatrix A1\n";
	cout << "-------------------------------------------------------------------\n";
	Print_matrix(A1);
	cout << "\n\nsolution = ( ";
	for (int i=0; i<N; i++) {
		cout << x[i] << " ";
	}
	cout << ")\n\ndet = " << Det(U) << endl;
}

int main() {
    ReadData();

    U = A;

    LU_d();
    Check();
    Solve();
    Inverse();

    WriteData();
    return 0;
}
