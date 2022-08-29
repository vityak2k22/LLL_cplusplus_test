#include "LLL.h"

// --------------------CORE OF PROGRAM--------------------

// LLL algorithm realization
void LLL(Matrix& X, double delta) {
    Matrix Y(X.size_m, X.size_n);
    Matrix coeff(X.size_m, X.size_m);
    Y.ortho_GS(X, coeff);
    SIZE k = 1;
    while (k < X.size_m) {
        for (SIZE j = k - 1; j >= 0; j--) {
            if (fabs(coeff.array[k][j]) > 0.5) {
				X.SubstractVectorWithMultiplier(round(coeff.array[k][j]), k, j);
                Y.ortho_GS(X, coeff);
            }
        }
        double a = (delta - coeff.array[k][k - 1] * coeff.array[k][k - 1]) * Y.Norm(k - 1);
        double b = Y.Norm(k);
        if (a <= b)
            k++;
        else {
            X.Swap(k, k - 1);
            Y.ortho_GS(X, coeff);
            k = (k - 1 > 1 ? k - 1 : 1);
        }
    }
}

// --------------------HELP METHODS--------------------

// Gram-Schmidt orthogonalization
void Matrix::ortho_GS(Matrix& X, Matrix& coeff) {
	AssignVector(X, 0);
    ARRAY B = new ELEMENT[size_m];
    B[0] = Norm(0);
    for (SIZE i = 1; i < size_m; i++) {
		AssignVector(X, i);
        for (SIZE j = 0; j < i; j++) {
            coeff.array[i][j] = ScalarProduct(*this, i, j) / B[j];
			SubstractVectorWithMultiplier(coeff.array[i][j], i, j);
        }
        B[i] = Norm(i);
    }
}

// Constructor for object creating
Matrix::Matrix(SIZE size_m, SIZE size_n) {
    this->size_m = size_m;
    this->size_n = size_n;
    array = new ARRAY [size_m];
    if (array == nullptr) {
        cout << "Memory 1 error\n";
        system("pause");
        exit(1);
    }
    for (SIZE i = 0; i < size_m; i++) {
        array[i] = new ELEMENT[size_n];
        if (array[i] == nullptr) {
            cout << i << "Memory 2 error\n";
            system("pause");
            exit(2);
        }
    }
    for (SIZE i = 0; i < size_m; i++) {
        for (SIZE j = 0; j < size_n; j++)
            array[i][j] = 0;
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& other) {
    size_m = other.size_m;
    size_n = other.size_n;
    array = new ARRAY [size_m];
    if (array == nullptr) {
        cout << "Memory 1 error\n";
        system("pause");
        exit(1);
    }
    for (SIZE i = 0; i < size_m; i++) {
        array[i] = new ELEMENT[size_n];
        if (array[i] == nullptr) {
            cout << i << "Memory 2 error\n";
            system("pause");
            exit(2);
        }
    }
    for (SIZE i = 0; i < size_m; i++) {
        for (SIZE j = 0; j < size_n; j++)
            array[i][j] = other.array[i][j];
    }
}

// Matrix input
void Matrix::write(istream &in) {
    for (SIZE i = 0; i < size_m; i++) {
        for (SIZE j = 0; j < size_n; j++)
            in >> array[i][j];
    }
}

// Matrix output
void Matrix::print(ostream &out) {
    for (SIZE i = 0; i < size_m; i++) {
        for (SIZE j = 0; j < size_n; j++)
            out << array[i][j] << " ";
        out << endl;
    }
}

// Matrix assignment
Matrix& Matrix::operator = (const Matrix& other) {
    if (this == &other)
        return *this;
    for (SIZE i = 0; i < size_m; i++) {
        for (SIZE j = 0; j < size_n; j++)
            array[i][j] = other.array[i][j];
    }
    return *this;
}

// Calculates scalar product vectors of matrix
double Matrix::ScalarProduct(Matrix Y, SIZE index1, SIZE index2) {
    double res = 0;
    for (SIZE j = 0; j < size_n; j++)
        res += array[index1][j] * Y.array[index2][j];
    return res;
}

// Calculates norm vectors of matrix
double Matrix::Norm(SIZE index) {
    double res = 0;
    for (SIZE j = 0; j < size_n; j++)
        res += array[index][j] * array[index][j];
    return res;
}

// Swaps vectors of matrix
void Matrix::Swap(SIZE index1, SIZE index2) {
    ARRAY temp = new ELEMENT[size_n];
    for (SIZE j = 0; j < size_n; j++) {
        temp[j] = array[index2][j];
        array[index2][j] = array[index1][j];
        array[index1][j] = temp[j];
    }
}

// Assigns basis vector values of another basis vector
void Matrix::AssignVector (Matrix& X, SIZE index) {
	for (SIZE j = 0; j < size_n; j++) {
		array[index][j] = X.array[index][j];
	}
}

// Assigns basis vector values according to the formula: [A = A - multiplier * B], where B is another basis vector
void Matrix::SubstractVectorWithMultiplier(double mult, SIZE index1, SIZE index2) {
	for (SIZE j = 0; j < size_n; j++) {
		array[index1][j] -= mult * array[index2][j];
	}
}

// Destructor
Matrix::~Matrix() {
    for (SIZE i = 0; i < size_m; i++) delete[] array[i];
    delete[] array;
}
