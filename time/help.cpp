#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double EPS = 1e-9;

// Функция перемножения двух матриц
vector<vector<double>> multiplyMatrices(const vector<vector<double>> &A, const vector<vector<double>> &B) {
    int n = A.size();
    int m = B[0].size();
    int p = A[0].size();
    vector<vector<double>> C(n, vector<double>(m, 0.0));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<double> sumVectors(const vector<double> &a, const vector<double> &b) {
    int n = a.size();
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}


// Функция умножения матрицы на вектор
vector<double> multiplyMatrixVector(const vector<vector<double>> &A, const vector<double> &b) {
    int n = A.size();
    vector<double> result(n, 0.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < A[i].size(); j++) {
            result[i] += A[i][j] * b[j];
        }
    }
    return result;
}

// Функция выполняет LUP-разложение
bool LUPDecomposition(vector<vector<double>> &A, vector<int> &P) {
    int n = A.size();
    for (int i = 0; i < n; i++) P[i] = i;
    
    for (int k = 0; k < n; k++) {
        double p = 0.0;
        int pivot = k;
        
        for (int i = k; i < n; i++) {
            if (fabs(A[i][k]) > p) {
                p = fabs(A[i][k]);
                pivot = i;
            }
        }
        
        if (p < EPS) return false; // Матрица вырождена
        
        swap(P[k], P[pivot]);
        swap(A[k], A[pivot]);
        
        for (int i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
    return true;
}

// Функция решает LUP * x = b
vector<double> LUPSolve(const vector<vector<double>> &A, const vector<int> &P, const vector<double> &b) {
    int n = A.size();
    vector<double> x(n), y(n);
    
    // Прямой ход (Ly = Pb)
    for (int i = 0; i < n; i++) {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= A[i][j] * y[j];
        }
    }
    
    // Обратный ход (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

vector<vector<double>> inverseMatrix(vector<vector<double>> A) {
    int n = A.size();
    vector<vector<double>> invA(n, vector<double>(n, 0.0));
    vector<int> P(n);
    
    if (!LUPDecomposition(A, P)) {
        throw runtime_error("Матрица вырождена, обратная матрица не существует.");
    }
    
    for (int i = 0; i < n; i++) {
        vector<double> e(n, 0.0);
        e[i] = 1.0;
        vector<double> col = LUPSolve(A, P, e);
        for (int j = 0; j < n; j++) {
            invA[j][i] = col[j];
        }
    }
    return invA;
}

// int main() {
//     int n;
//     cout << "Введите размерность системы: ";
//     cin >> n;
    
//     vector<vector<double>> A(n, vector<double>(n));
//     vector<double> b(n);
//     vector<int> P(n);
    
//     cout << "Введите коэффициенты матрицы A:" << endl;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             cin >> A[i][j];
//         }
//     }
    
//     cout << "Введите вектор b:" << endl;
//     for (int i = 0; i < n; i++) {
//         cin >> b[i];
//     }
    
//     if (!LUPDecomposition(A, P)) {
//         cout << "Матрица вырождена, решение невозможно." << endl;
//         return 1;
//     }
    
//     vector<double> x = LUPSolve(A, P, b);
    
//     cout << "Решение системы:" << endl;
//     for (double xi : x) {
//         cout << xi << " ";
//     }
//     cout << endl;
    
//     return 0;
// }
