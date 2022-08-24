#include "LLL.h"

int main () {
    ifstream in("input.txt");
    ofstream out("output.txt");
    
    out << "LLL algorithm (my code describing <double>)\n\n";
    long long m, n;
    in >> m >> n;
    Matrix X(m, n), Y(m, n);
    out << "A matrix " << m << "x" << n << endl;
    X.write(in);
    out << "Input matrix:\n";
    X.print(out);
    LLL(Y);
    out << "\nReduced matrix:\n";
    Y.print(out);
    in.close();
    out.close();
	system("pause");
	return 0;
}
