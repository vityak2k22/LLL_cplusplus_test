#ifndef LLL_H
#define LLL_H


#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>

using namespace std;

class Matrix;
void LLL(Matrix& X, double delta = 0.75);

class Matrix {
    double** array;
    long long size_m, size_n;
public:
	Matrix(long long size_m, long long size_n);
	Matrix(const Matrix& other);
    void write(istream &in);
	void print(ostream &out);
	Matrix& operator = (const Matrix& other);
	double ScalarProduct(Matrix Y, long long index1, long long index2);
	double Norm(long long index);
	void Swap(long long index1, long long index2);
	Matrix ortho_GS(Matrix& coeff);
	friend void LLL(Matrix& X, double delta);
	~Matrix();
};


#endif
