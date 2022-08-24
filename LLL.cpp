#include "LLL.h"



Matrix::Matrix(long long size_m, long long size_n) {
    this->size_m = size_m;
    this->size_n = size_n;
    array = new double* [size_m];
    if (array == nullptr) {
        cout << "Memory 1 error\n";
        system("pause");
        return;
    }
    for (long long i = 0; i < size_m; i++) {
        array[i] = new double[size_n];
        if (array[i] == nullptr) {
            cout << i << "Memory 2 error\n";
            system("pause");
            return;
        }
    }
    for (long long i = 0; i < size_m; i++) {
        for (long long j = 0; j < size_n; j++)
            array[i][j] = 0;
    }
}

Matrix::Matrix(const Matrix& other) {
    size_m = other.size_m;
    size_n = other.size_n;
    array = new double* [size_m];
    if (array == nullptr) {
        cout << "Memory 1 error\n";
        system("pause");
        return;
    }
    for (long long i = 0; i < size_m; i++) {
        array[i] = new double[size_n];
        if (array[i] == nullptr) {
            cout << i << "Memory 2 error\n";
            system("pause");
            return;
        }
    }
    for (long long i = 0; i < size_m; i++) {
        for (long long j = 0; j < size_n; j++)
            array[i][j] = other.array[i][j];
    }
}

void Matrix::write(istream &in) {
    for (long long i = 0; i < size_m; i++) {
        for (long long j = 0; j < size_n; j++)
            in >> array[i][j];
    }
}

void Matrix::print(ostream &out) {
    out << "[";
    for (long long i = 0; i < size_m; i++) {
        out << "[";
        for (long long j = 0; j < size_n; j++) {
            if (j == size_n - 1)
                out << array[i][j];
            else
                out << array[i][j] << ", ";
        }
        if (i == size_m - 1)
            out << "]";
        else
            out << "], ";
    }
    out << "]";
}

Matrix& Matrix::operator = (const Matrix& other) {
    if (this == &other)
        return *this;
    for (long long i = 0; i < size_m; i++) {
        for (long long j = 0; j < size_n; j++)
            array[i][j] = other.array[i][j];
    }
    return *this;
}

double Matrix::ScalarProduct(Matrix Y, long long index1, long long index2) {
    double res = 0;
    for (long long j = 0; j < size_n; j++)
        res += array[index1][j] * Y.array[index2][j];
    return res;
}

double Matrix::Norm(long long index) {
    double res = 0;
    for (long long j = 0; j < size_n; j++)
        res += array[index][j] * array[index][j];
    return res;
}

void Matrix::Swap(long long index1, long long index2) {
    double* temp = new double[size_n];
    for (long long j = 0; j < size_n; j++) {
        temp[j] = array[index2][j];
        array[index2][j] = array[index1][j];
        array[index1][j] = temp[j];
    }
}

Matrix Matrix::ortho_GS(Matrix& coeff) {
    Matrix Y(size_m, size_n);
    for (long long x = 0; x < size_n; x++)
        Y.array[0][x] = array[0][x];
    double* B = new double[size_m];
    B[0] = Y.Norm(0);
    for (long long i = 1; i < size_m; i++) {
        for (long long x = 0; x < size_n; x++)
            Y.array[i][x] = array[i][x];
        for (long long j = 0; j < i; j++) {
            coeff.array[i][j] = ScalarProduct(Y, i, j) / B[j];
            for (long long x = 0; x < size_n; x++)
                Y.array[i][x] -= coeff.array[i][j] * Y.array[j][x];
        }
        B[i] = Y.Norm(i);
    }
    return Y;
}

void LLL(Matrix& X, double delta) {
    Matrix Y(X.size_m, X.size_n);
    Matrix coeff(X.size_m, X.size_m);
    Y = X.ortho_GS(coeff);
    long long k = 1;
    while (k < X.size_m) {
        for (long long j = k - 1; j >= 0; j--) {
            if (fabs(coeff.array[k][j]) > 0.5) {
                for (long long x = 0; x < X.size_n; x++)
                    X.array[k][x] -= round(coeff.array[k][j]) * X.array[j][x];
                Y = X.ortho_GS(coeff);
            }
        }
        double a = (delta - coeff.array[k][k - 1] * coeff.array[k][k - 1]) * Y.Norm(k - 1);
        double b = Y.Norm(k);
        if (a <= b)
            k++;
        else {
            X.Swap(k, k - 1);
            Y = X.ortho_GS(coeff);
            k = (k - 1 > 1 ? k - 1 : 1);
        }
    }
}

Matrix::~Matrix() {
    for (long long i = 0; i < size_m; i++) delete[] array[i];
    delete[] array;
}
