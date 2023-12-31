#ifndef VECTOR_ARITHMETICS_HPP
#define VECTOR_ARITHMETICS_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;


// OPERATORS.
vector<double> operator+(const vector<double>& u, const vector<double>& v){
    /* Addition operator between two vectors. */
    assert(u.size() == v.size());
    vector<double> res(u.size(), 0.);

    for (size_t i = 0; i < u.size(); i++){
        res[i] = u[i] + v[i];
    }

    return res;
}

vector<double> operator-(const vector<double>& u, const vector<double>& v){
    /* Substraction operator between two vectors. */
    assert(u.size() == v.size());
    vector<double> res(u.size(), 0.);

    for (size_t i = 0; i < u.size(); i++){
        res[i] = u[i] - v[i];
    }

    return res;
}

void operator+=(vector<double>& u, const vector<double>& v){
    /* Self-Add operator between two vectors. */
    assert(u.size() == v.size());
    for (size_t i = 0; i < u.size(); i++){
        u[i] += v[i];
    }
}

vector<double> operator*(const vector<double>& u, const double& l){
    /* Multiplication operator between a vector and a scalar. */
    vector<double> res(u.size(), 0.);
    for (size_t i = 0; i < u.size(); i++){
        res[i] = u[i] * l;
    }

    return res;
}

vector<double> operator*(const double& l, const vector<double>& u){
    /* Multiplication operator between a vector and a scalar. */
    vector<double> res(u.size(), 0.);
    for (size_t i = 0; i < u.size(); i++){
        res[i] = u[i] * l;
    }

    return res;
}

vector<double> operator/(const vector<double>& u, const double& l){
    /* Division operator between a vector and a scalar. */
    vector<double> res(u.size(), 0.);
    for (size_t i = 0; i < u.size(); i++){
        res[i] = u[i] / l;
    }

    return res;
}

ostream& operator<<(ostream& o, const vector<double>& u){
    /* Stream operator. */
    o << "[";
    for (size_t i = 0; i < u.size(); i++){
        o << u[i];
        if (i != u.size() - 1){
            o << "\t";
        }
    }
    o << "]" << endl;
    
    return o;
}


// FUNCTIONS.
double scalarProduct(const vector<double>& u, const vector<double>& v){
    /* Compute the scalar product between two vectors. */
    assert(u.size() == v.size());
    double res = 0;

    for (size_t i = 0; i < u.size(); i++){
        res += u[i] * v[i];
    }

    return res;
}

double Norm2(const vector<double>& u){
    /* Compute the Norm 2 of a vector. */
    double res = 0;
    for (size_t i = 0; i < u.size(); i++){
        res += u[i] * u[i];
    }
    return sqrt(res);

}



#endif