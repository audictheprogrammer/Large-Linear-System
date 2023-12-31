#ifndef MAP_MATRIX_HPP
#define MAP_MATRIX_HPP

using namespace std;
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

template <typename type_t>
class MapMatrix{
    public:
	int nr;
	int nc;
	typedef tuple<int, int> NxN;
	map<NxN, type_t> data;

	MapMatrix(const int& nr=0, const int& nc=0){
		/* Constructor. */
		this->nr = nr;
		this->nc = nc;
	}

	MapMatrix(const MapMatrix& M){
		/* Constructor. */
		this->nr = M.nr;
		this->nc = M.nc;
		this->data = M.data;
	}

	MapMatrix operator =(const MapMatrix& M){
		/* Equal operator. */
		return MapMatrix(M);
	}

	void insert(const int& i, const int& j, const type_t& v){
		// Checking if j,k are correct indices.
		assert(0<=i && i<nr && 0<=j && j<nc);

		// Inserting.
		NxN key = {i, j};
		this->data[key] = v;
	}

	double& operator()(const int& i, const int& j){
		NxN key = {i, j};
		return data[key];
	}

	void operator +=(MapMatrix& M){
		/* Self-adds with a MapMatrix. */
		assert(this->nc == M.nc && this->nr == M.nr);

		for (auto it = M.data.begin(); it != M.data.end(); it++){
			NxN key = it->first;
			this->data[key] += M.data[key];
		}
	}

	void operator *=(MapMatrix& M2){
		/* Self-multiplies with a MapMatrix.
        Can probably be done in nlogn.
        */
		MapMatrix& M1 = *this;
		assert(M1.nc == M2.nr);

		MapMatrix res(M1.nr, M2.nc);
		for (auto it = M1.data.begin(); it != M1.data.end(); it++){
			for (auto it2 = M2.data.begin(); it2 != M2.data.end(); it2++){
				int i1 = get<0> (it->first);
				int j1 = get<1> (it->first);
				int i2 = get<0> (it2->first);
				int j2 = get<1> (it2->first);
				if (j1 == i2){
					NxN key = {i1, j2};
					res.data[key] += it->second * it2->second;
				}
			}
		}
        *this = MapMatrix(res);

	}

	void operator *=(type_t l){
		/* Self-multiplies with a type_t. */
		for (auto it = this->data.begin(); it != this->data.end(); it++){
			// this->data[it->first] *= l;
			it->second *= l;  // Cela fonctionne;
		}
	}

	MapMatrix operator +(MapMatrix& M2){
		/* Adds two MapMatrix. */
		MapMatrix& M1 = *this;
		assert(M1.nc == M2.nc && M1.nr == M2.nr);

        MapMatrix res(M1.nr, M1.nc);
        for (auto it = M1.data.begin(); it != M1.data.end(); it++){
			NxN key = it->first;
			res.data[key] += M1.data[key];
		}

        for (auto it = M2.data.begin(); it != M2.data.end(); it++){
            NxN key = it->first;
            res.data[key] += M2.data[key];
        }

		return res;
	}

	MapMatrix operator *(MapMatrix& M2){
		/* Multiplies two MapMatrix. */
		MapMatrix& M1 = *this;
		assert(M1.nc == M2.nr);

		MapMatrix res(M1.nr, M2.nc);
		for (auto it = M1.data.begin(); it != M1.data.end(); it++){
			for (auto it2 = M2.data.begin(); it2 != M2.data.end(); it2++){
				int i1 = get<0> (it->first);
				int j1 = get<1> (it->first);
				int i2 = get<0> (it2->first);
				int j2 = get<1> (it2->first);
				if (j1 == i2){
					// NxN key = {i1, j2};
					// res.data[key] += it->second * it2->second;
					res.insert(i1, j2, it->second* it2->second);
				}
			}
		}
		return res;
	}

	vector<type_t> operator *(vector<type_t> V){
		/* Multiplies a MapMatrix with a vector. */
		assert(((size_t) this->nc) == V.size());
		vector<type_t> res(this->nr);
		for (auto it = this->data.begin(); it != this->data.end(); it++){
			int i = get<0> (it->first);
			int j = get<1> (it->first);
			res[i] +=  this->data[it->first] * V[j];
		}
		return res;
	}


	friend const int& NbRow(const MapMatrix& M){
		return M.nr;
	}
	friend const int& NbCol(const MapMatrix& M){
		return M.nc;
	}


};


template<typename type_t>
ostream& operator <<(ostream& o, const MapMatrix<type_t>& M){
	/* Output stream operator. */
	for (auto it = M.data.begin(); it != M.data.end(); it++){
		tuple<int, int> key = it->first;
		int i = get<0> (key);
		int j = get<1> (key);
		o << "(" << i << "," << j << "): \t" << it->second << "\n";
	}

	return o;
}

template<typename type_t>
MapMatrix<type_t> operator *(type_t l, MapMatrix<type_t>& M){
	/* Operator: Multiplies a type_t (left) with a MapMatrix (right). */
	MapMatrix<type_t> res(M.nr, M.nc);
	for (auto it = M.data.begin(); it != M.data.end(); it++){
		res.data[it->first] = l * it->second;
	}
	return res;
}

MapMatrix<double> LoadMapMatrix(const string& filename){
	/* Reads a file and construct the MapMatrix. */
    fstream f;
    f.open(filename);

	if (f.is_open() == false){
		cout << "File is not open !" << endl;
		return MapMatrix<double> (0, 0);
	}
	string line;
	istringstream iss;

	if (!getline(f, line)){
		cout << "Line is not read !" << endl;
		return MapMatrix<double> (0, 0);
	}

	istringstream ss(line);
	int nr;
	int nc;
	ss >> nr >> nc;
	cout << "NR = " << nr << endl << "NC = " << nc << endl;
	MapMatrix<double> res(nr, nc);

	int row, col;
	double val;
	while (getline(f, line)) {
		istringstream ss(line); // input string -> stream;
		ss >> row >> col >> val;
		res.insert(row, col, val);
		iss.clear();
	}
	f.close();
	return res;
}


void Write(const string& filename, const MapMatrix<double>& m){
	/* Writes an matrix into a file. */
	ofstream f;
	f.open(filename);

	if (!f.is_open()){
		cout << "Error f.open: NOT OPEN !" << endl;
	}
	f << m.nr << " " << m.nc << endl;
	for (auto it = m.data.begin(); it != m.data.end(); it++){
		int i = get<0> (it->first);
		int j = get<1> (it->first);
		double value = it->second;

		f << i << " " << j << " " << value << endl;
		cout << i << " " << j << " " << value << endl;
	}

	f.close();

}




#endif
