#ifndef DENSE_MATRIX_HPP
#define DENSE_MATRIX_HPP

#include <cassert>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>


using namespace std;

class DenseMatrix{
public:
    int nr;  // number of rows;
    int nc;  // number of columns;
    vector<double> data;  // Vector storing the entries row-wise;

DenseMatrix(const int& nr=0, const int& nc=0)
:nr(nr), nc(nc), data(nc*nr, 0.) {}


DenseMatrix(const DenseMatrix& M){
    this->nr = M.nr;
    this->nc = M.nr;
    this->data = M.data;
}

DenseMatrix operator =(DenseMatrix M){
    return DenseMatrix(M);
}

double& operator()(const int& i , const int& j){
    /* Internal operator for evaluation. */
    return this->data[i * nc + j];
}

double operator()(const int& i , const int& j) const{
    /* External operator for evaluation. */
    return this->data[i * nc + j];
}


vector<double> row(const int& i){
    /* Returns the row i of the DenseMatrix. */
    vector<double> res(this->nr);
    for (size_t j = 0; j < res.size(); j++){
        res[j] = (*this)(i, j);
    }
    return res;
}

vector<double> col(const int& j){
    /* Returns the column j of the DenseMatrix. */
    vector<double> res(this->nc);
    for (size_t i = 0; i < res.size(); i++){
        res[i] = (*this)(i, j);
    }
    return res;
}

DenseMatrix operator +(DenseMatrix& M2){
	/* Sums a DenseMatrix with a DenseMatrix. */
	DenseMatrix& M1 = *this;
	assert(M1.nr == M2.nr && M1.nc == M2.nc);
	DenseMatrix res(M1.nc, M1.nr);
	for (int i = 0; i < res.nc * res.nr; i++){
		res.data[i] = M1.data[i] + M2.data[i];
	}

	return res;
}

void operator +=(DenseMatrix& M2){
	/* Sums a DenseMatrix with a DenseMatrix and modifies the first matrix. */
	DenseMatrix& M1 = *this;
	assert(M1.nr == M2.nr && M1.nc == M2.nc);
	for (int i = 0; i < M1.nc * M1.nr; i++){
		M1.data[i] += M2.data[i];
	}

}

DenseMatrix operator*(DenseMatrix& M2){
	/* Multiplies a DenseMatrix with a DenseMatrix. */
	DenseMatrix& M1 = *this;
	assert(M1.nc == M2.nr);
	DenseMatrix res(M1.nr, M2.nc);

	for (int i = 0; i < M1.nr; i++){
		for (int j = 0; j < M2.nc; j++){
			for (int k = 0 ; k < M1.nc; k++){
				res.data[i*res.nc + j] += M1.data[i*M1.nc + k] * M2.data[k*M2.nc + j];
			}
		}
	}

	return res;
}


vector<double> operator*(vector<double> V){
	/* Multiplies a DenseMatrix by a vector. */
	assert(this->nc == (int) V.size());
	vector<double> res(this->nr);
	for (int i = 0; i < this->nr; i++){
		for (int j = 0; j < this->nc; j++){
			res[i] += this->data[i*this->nc + j] * V[j];
		}
	}

	return res;
}

void operator*=(DenseMatrix& M2){
	/* Multiplies a DenseMatrix with a DenseMatrix and modifies the first matrix.*/
	*this = *this*M2;
}


void operator*=(const double l){
	/* Multiplies a DenseMatrix with a double. */
	for (int i = 0; i < this->nc*this->nr; i++){
		this->data[i] *= l;
	}
}


friend int NbRow(const DenseMatrix& M){
	return M.nr;
}

friend int NbCol(const DenseMatrix& M){
	return M.nc;
}

};


DenseMatrix operator*(const double& l, DenseMatrix& M){
	/* Multiplies a double(left) with a matrix (right). */
	DenseMatrix res(M.nr, M.nc);
	for (int i = 0; i < M.nr * M.nc; i++){
		res.data[i] = l * M.data[i];
	}
	return res;
}

ostream& operator<<(ostream& o, const DenseMatrix& M){
    /* Output stream operator. */
    for (int i = 0; i < M.nr * M.nc; i++){
    	if ((i+1)%M.nc == 0){
    		o << std::right << M.data[i] << "\n";
    	} else {
    		o << std::right << M.data[i] << "\t";
    	}
    }

    return o;
}


void LUFactorize(DenseMatrix& M){
    /* Factorize the initial matrix into a LU form.
    Does not use any permutation: M = LU. */
    int N = M.nr;
    for (int i = 0; i < N-1; i++){
        double pivot = M(i, i);
        assert(pivot != 0);

        for (int j = i+1; j < N; j++){
            M(j, i) /= pivot;  // Lower part
            for (int k = i+1; k < N; k++){
                M(j, k) -= M(i, k) * M(j, i); // Upper part
            }
        }
    }

}


void LUFactorize(DenseMatrix& M, vector<int>& P){
    /* Factorize the initial matrix into a LU form.
    Using a pre-permutation on M: P*M = LU.
    We consider the permutation is set to Identity i.e. P = [0, .., N-1]
    */
    size_t N = M.nr;
    assert(P.size() == N);
    
    for (size_t i = 0; i < N-1; i++){

        // Finding the pivot.
        int index = i;
        double pivot = M(P[i], i);
        for (size_t j = i+1; j < N; j++){
            if (std::abs(pivot) < std::abs(M(P[j], i))) {
                index = j;
                pivot = M(P[j], i);
            }
        }
        std::swap(P[i], P[index]);

        // Performing the gaussian elimination.
        for (size_t j = i+1; j < N; j++){
            M(P[j], i) /= pivot;  // Lower part
            for (size_t k = i+1; k < N; k++){
                M(P[j], k) -= M(P[i], k) * M(P[j], i); // Upper part
            }
        }

    }

}

DenseMatrix LoadDenseMatrix(const string& filename){
    /* Load from a file a DenseMatrix. */
    ifstream file(filename);
    assert(file.is_open() == 1);

    string str;
    // Line 1: #SIZE.
    file >> str;
    assert("#SIZE" == str);

    // Line 2 & 3: nr and nc.
    int nr, nc;
    file >> nr;
    file >> nc;
    cout << "nr = " << nr << endl;
    cout << "nc = " << nc << endl;

    DenseMatrix res(nr, nc);

    // Line 4: #DATA.
    file >> str;
    assert("#DATA" == str);

    // Next lines: Matrix coefficients.
    int i = 0, j = 0;
    while(file){
        file >> res(i, j);

        j++;
        if (j == nc){
            j = 0;
            i++;
        }
    }

    return res;
}



#endif
