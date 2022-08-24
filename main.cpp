#include "LLL.h"

int main () {
	out << "LLL algorithm (my code describing <double>)\n\n";
    long long m, n;
    in >> m >> n;
    Matrix X(m, n), Y(m, n);
    out << "A matrix " << m << "x" << n << endl;
    X.write();
    out << "Input matrix:\n";
    X.print();
    LLL(Y);
    out << "\nReduced matrix:\n";
    Y.print();
    in.close();
    out.close();
	system("pause");
	return 0;
}