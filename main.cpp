#include <iostream>
#include <vector>
#include "MapMatrix.hpp"
#include "DenseMatrix.hpp"
#include "FredholmMatrix.hpp"
#include "VectorArithmetics.hpp"

using namespace std;
int main(){

    /* Testing the minimal residual linear solving function.
        [1 0 0]             [1]
    A = [0 1 0] + [1 1 1] * [1]
        [0 0 1]             [1]
    */
    vector<double> u1 = {1, 1, 1};
    vector<double> v1 = {1, 1, 1};

    FredholmMatrix A(3);
    A.insert(0, 0, 1);
    A.insert(1, 1, 1);
    A.insert(2, 2, 1);
    A.insert(u1, v1);

    cout << A << endl;
    vector<double> b = {1, 2, 3};
    vector<double> x = MinResSolve(A, b);

    cout << x << endl;


    // std::cout.precision(5);
    // PlotGraph(10);

    // std::cout.precision(5);
    // PlotGraph("matrix.txt");


    return 0;
}
