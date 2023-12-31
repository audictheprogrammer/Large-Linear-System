#ifndef FREDHOLM_MATRIX_HPP
#define FREDHOLM_MATRIX_HPP


using namespace std;
#include <map>
#include <vector>
#include "MapMatrix.hpp"
#include "VectorArithmetics.hpp"
#include <cassert>

class FredholmMatrix{
    public:
    int n; // Size of the square matrix
    MapMatrix<double> D; // A sparse matrix modelling the contribution of D
    vector<vector<double>> lru; // Representing the vectors （u_j）
    vector<vector<double>> lrv; // Representing the vectors （v_j）

    // CONSTRUCTORS.
    FredholmMatrix(const int& n0 = 0): n(n0), D(n0, n0), lru(), lrv() {}

    FredholmMatrix(const FredholmMatrix& M): n(M.n), D(M.D), lru(M.lru), lrv(M.lrv) {}

    // OPERATORS.
    FredholmMatrix& operator=(const FredholmMatrix& M){
        /* Equal operator. */
        *this = FredholmMatrix(M);
        return *this;
    }

    vector<double> operator*(const vector<double>& u){
        /* Matrix-Vector product operator.
        Complexity: O(r * n).
        */
        assert(u.size() == ((size_t) n));
        vector<double> res(n);

        // Performing the D*u product.
        res = D * u;

        // Performing the sum: Sum[u_j * v_j^T * u] = Sum[u_j * scalar].
        size_t r = lru.size();
        for (size_t j = 0; j < r; j++){
            double scalar = scalarProduct(lrv[j], u);
            for (int i = 0; i < n; i++){
                res[i] += lru[j][i] * scalar;
            }

        }

        return res;
    }

    double operator()(const int& i, const int& j){
        /* Evaluation operator. Returns a copy. */
        double res = 0;
        res += D(i, j);

        size_t r = lru.size();
        for (size_t k = 0; k < r; k++){
            res += lru[k][i] * lrv[k][j];
        }

        return res;
    }

    // FUNCTIONS.
    void insert(const int& i, const int& j, const double& val){
        /* Insert an element into the FredholmMatrix. */
        assert(0 <= i && i < this->n);
        assert(0 <= j && j < this->n);
        this->D(i, j) += val;
    }

    void insert(const vector<double>& u, const vector<double>& v){
        /* Insert u and v into the FredholmMatrix. */
        this->lru.push_back(u);
        this->lrv.push_back(v);
    }

    friend vector<double> MinResSolve(const FredholmMatrix& A, const vector<double> b);
};

ostream& operator<<(ostream& o, FredholmMatrix& M){
    /* Stream operator. */
    o << "[";
    for (int i = 0; i < M.n; i++){
        for (int j = 0; j < M.n; j++){
            double v = M(i, j);
            o << v;
            if (j == M.n-1){
                o << endl;
            } else {
                o << "\t";
            }
        }
    }
    o << "]" << endl;

    return o;
}

vector<double> MinResSolve(FredholmMatrix& A, vector<double>& b){
    /* Solve the linear system Ax = b using the minimal residual method.
    r = b - Ax
    a = r^T A r / Norm(Ar, 2)^2
    x = x + ar
    */
    assert(b.size() == ((size_t) A.n));
    vector<double> r(A.n, 0.);
    double a = 0;
    vector<double> x(A.n, 0.);

    int nb = 0;
    while (nb < 10){
        r = b - A*x;
        vector<double> tmp = A*r;
        a = scalarProduct(r, tmp) / (Norm2(tmp) * Norm2(tmp));
        x = x + a*r;
    }

    return x;
}



#endif

