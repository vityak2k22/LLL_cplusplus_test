#ifndef LLL_H
#define LLL_H


#include <iostream>
#include <fstream>
#include <cmath>

typedef long long SIZE;
typedef double** MATRIX;
typedef double* ARRAY;
typedef double ELEMENT;

using namespace std;

// Matrix class is created for work with lattice basis (prototype)
class Matrix;
// LLL algorithm realization (prototype of function)
void LLL(Matrix& X, double delta = 0.75);

// Matrix class is created for work with lattice basis
class Matrix {
    MATRIX array;																// Array of basis vectors. Presented as a matrix of coordinates
    SIZE size_m, size_n; 														// Dimensions of matrix
public:
	Matrix(SIZE size_m, SIZE size_n);											// Constructor for object creating
	Matrix(const Matrix& other);												// Copy constructor
    void write(istream &in);													// Matrix input
	void print(ostream &out);													// Matrix output
	friend void LLL(Matrix& X, double delta);									// LLL algorithm realization - CORE OF PROGRAM
	~Matrix();
private:
	Matrix ortho_GS(Matrix& coeff);												// Gram-Schmidt orthogonalization
	Matrix& operator = (const Matrix& other);									// Matrix assignment
	void AssignVector (Matrix& X, SIZE index);									// Assigns basis vector values of another basis vector
	void SubstractVectorWithMultiplier(double mult, SIZE index1, SIZE index2);	// Assigns basis vector values according to the formula: [A = A - multiplier * B], where B is another basis vector
	double ScalarProduct(Matrix Y, SIZE index1, SIZE index2);					// Calculates scalar product vectors of matrix
	double Norm(SIZE index);													// Calculates norm vectors of matrix
	void Swap(SIZE index1, SIZE index2);										// Swaps vectors of matrix
};


#endif
