#include <iostream>

using namespace std;

#define N 2

class Complex {
private:
    double real;
    double imag;

public:
    Complex(double r = 0.0, double i = 0.0) : real(r), imag(i) {}

    double getReal() const { return real; }
    double getImag() const { return imag; }

    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator*(const Complex& other) const {
        return Complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    Complex operator/(const Complex& other) const {
        double denominator = other.real * other.real + other.imag * other.imag;
        return Complex((real * other.real + imag * other.imag) / denominator, (imag * other.real - real * other.imag) / denominator);
    }

    void print() const {
        cout << real << " + " << imag << "i";
    }
};

typedef Complex Matrix[N][N];
typedef Complex Vector[N];

// Function to perform LU decomposition
void luDecomposition(const Matrix& A, Matrix& L, Matrix& U) {
    for (int i = 0; i < N; i++) {
        // Upper Triangular Matrix U
        for (int k = i; k < N; k++) {
            Complex sum(0, 0);
            for (int j = 0; j < i; j++) {
                sum = sum + L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - sum;
        }

        // Lower Triangular Matrix L
        for (int k = i; k < N; k++) {
            if (i == k)
                L[i][i] = Complex(1, 0);  // Diagonal as 1
            else {
                Complex sum(0, 0);
                for (int j = 0; j < i; j++) {
                    sum = sum + L[k][j] * U[j][i];
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Function to solve system of equations using forward and backward substitution
void solveLU(const Matrix& L, const Matrix& U, const Vector& b, Vector& x) {
    Complex y[N];

    // Forward substitution to solve L * y = b
    for (int i = 0; i < N; i++) {
        Complex sum(0, 0);
        for (int j = 0; j < i; j++) {
            sum = sum + L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    // Backward substitution to solve U * x = y
    for (int i = N - 1; i >= 0; i--) {
        Complex sum(0, 0);
        for (int j = i + 1; j < N; j++) {
            sum = sum + U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
}

// Function to calculate the inverse of a matrix using LU decomposition
void inverseMatrix(const Matrix& A, Matrix& inverse) {
    Matrix L, U;
    luDecomposition(A, L, U);

    for (int i = 0; i < N; i++) {
        Vector e = {Complex(0, 0)};
        e[i] = Complex(1, 0);
        Vector column;
        solveLU(L, U, e, column);
        for (int j = 0; j < N; j++) {
            inverse[j][i] = column[j];
        }
    }
}

// Function to print a matrix
void printMatrix(const Matrix& matrix) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j].print();
            cout << " ";
        }
        cout << endl;
    }
}

int main() {
    // Example of complex matrix
    Matrix A = {
        {Complex(2, 1), Complex(3, -1)},
        {Complex(1, 2), Complex(4, 0)}
    };

    cout << "Original Matrix A:" << endl;
    printMatrix(A);

    Matrix A_inv;
    inverseMatrix(A, A_inv);

    cout << "\nInverse Matrix A^-1:" << endl;
    printMatrix(A_inv);

    return 0;
}
