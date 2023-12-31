#ifndef FREDHOLM_MATRIX_HPP
#define FREDHOLM_MATRIX_HPP


using namespace std;
#include <map>
#include <vector>
#include "MapMatrix.hpp"
#include "VectorArithmetics.hpp"
#include <cassert>
#include <fstream>


#define EPS 1.0E-10
#define NB_ITER 100
#define MAX_R 15
#define MAX_R_2 100


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

};

ostream& operator<<(ostream& o, FredholmMatrix& M){
    /* Stream operator. */
    o << "[";
    for (int i = 0; i < M.n; i++){
        for (int j = 0; j < M.n; j++){
            double v = M(i, j);
            o << v;
            if (j != M.n-1){
                o << "\t";
            } else if (i != M.n - 1){
                o << endl;
            }
        }
    }
    o << "]"<< endl;

    return o;
}

vector<double> MinResSolve(FredholmMatrix& A, vector<double>& b){
    /* Solves the linear system Ax = b using the minimal residual method.
    r = b - Ax
    a = r^T A r / Norm(Ar, 2)^2
    x = x + ar
    */
    assert(b.size() == ((size_t) A.n));
    vector<double> r(A.n, 0.);
    double a = 0;
    vector<double> x(A.n, 0.);

    int nb = 0;
    while (nb < NB_ITER){
        r = b - A*x;
        vector<double> tmp = A*r;
        a = scalarProduct(r, tmp) / (Norm2(tmp) * Norm2(tmp));
        x = x + a*r;
        nb++;
        if (Norm2(a*r) < EPS){
            cout << "MinResSolve completed with nb iter: " << nb << endl;
            return x;
        }
    }

    cout << "MinResSolve completed with nb iter: " << nb << endl;
    return x;
}

FredholmMatrix CrossApproximation(const DenseMatrix& B, const int& r){
    /* Computes a low-rank approximation of B. 
    */
    assert(B.nc == B.nr);

    DenseMatrix R = B;
    FredholmMatrix Bt(B.nr); // Representing B_tilda

    int j0 = 0, k0 = 0;
    double mx = std::abs(R(j0, k0));
    for (int p = 0; p < r; p++){
        // Set j0, k0 to argmax.
        for (int j = 0; j < R.nr; j++){
            for (int k = 0; k < R.nc; k++){
                double l = std::abs(R(j, k));
                // cout << "L = " << l << endl;
                if (mx < l){
                    mx = l;
                    j0 = j;
                    k0 = k;
                }
            }
        }

        // Get up and vp.
        vector<double> up = R.col(k0) / R(j0, k0);
        vector<double> vp = R.row(j0);

        // Move u_p and v_p from R to Bt.
        for (int i = 0; i < R.nr; i++){
            for (int j = 0; j < R.nc; j++){
                R(i, j) -= up[i] * vp[j]; // up[i] or vp[i]
            }
        }

        Bt.insert(up, vp);

        mx = 0;
    }


    return Bt;
}


double FrobeniusNorm(const DenseMatrix& M){
    /* Computes the Frobenius norm of a DenseMatrix. */
    double res = 0.;

    for (int i = 0; i < M.nr; i++){
        for (int j = 0; j < M.nc; j++){
            res += M(i, j) * M(i, j);
        }
    }
    return std::sqrt(res);
}

void Write(const string& filename, const vector<double>& X, const vector<double>& Y){
    /* Writes data into a file. */
    ofstream f;
    f.open(filename);

    if (!f.is_open()){
        cout << "Error f.open: NOT OPEN !" << endl;
    }
        for (size_t i = 0; i < X.size(); i++){
            f << X[i] << " " << Y[i] << endl;
            cout << X[i] << " " << Y[i] << endl;
        }

        f.close();
    }

void PlotGraph(const size_t& n){
    /* Plot the function log(r) -> log(FrobNorm(B - Br)).
    B = (exp[-(j-k)^2 / n^2].
    */

    // The X_axis and Y_axis for plotting.
    vector<double> X(MAX_R);
    vector<double> Y(MAX_R);

    // Constructing B.
    DenseMatrix B(n, n);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < n; j++){
            double t = -std::pow((double)i-j, 2) / std::pow(n, 2);
            double v = std::exp(t);
            B(i, j) = v;
        }
    }

    // Constructing Br.
    for (int r = 1; r <= MAX_R; r++){
        cout << "r = " << r << endl;
        FredholmMatrix Br = CrossApproximation(B, r);

        // Computing the norm of B-Br.
        DenseMatrix Diff(n, n);
        for (size_t i = 0; i < n; i++){
            for (size_t j = 0; j < n; j++){
                Diff(i, j) = B(i, j) - Br(i, j);
            }
        }
        X[r-1] = r;
        Y[r-1] = FrobeniusNorm(Diff);
    }

    Write("Graph1.txt", X, Y);

}

void PlotGraph(const string& filename){
    /* Plot the function log(r) -> log(FrobNorm(B - Br)).
    B loaded from filename.
    */
    
    // The X_axis and Y_axis for plotting.
    vector<double> X(MAX_R_2);
    vector<double> Y(MAX_R_2);

    // Loading B.
    DenseMatrix B = LoadDenseMatrix(filename);

    // Constructing Br.
    for (int r = 1; r <= MAX_R_2; r++){
        cout << "r = " << r << endl;
        FredholmMatrix Br = CrossApproximation(B, r);

        // Computing the norm of B-Br.
        DenseMatrix Diff(B.nr, B.nc);
        for (int i = 0; i < B.nr; i++){
            for (int j = 0; j < B.nc; j++){
                Diff(i, j) = B(i, j) - Br(i, j);
            }
        }
        X[r-1] = r;
        Y[r-1] = FrobeniusNorm(Diff);
    }

    Write("Graph2.txt", X, Y);
}

#endif

