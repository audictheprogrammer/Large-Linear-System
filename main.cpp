#include <iostream>
#include <vector>
#include "MapMatrix.hpp"
#include "DenseMatrix.hpp"
#include "FredholmMatrix.hpp"
#include "VectorArithmetics.hpp"
#include <cstdlib>

using namespace std;


void input(int options[], int argc, char* argv[]){
    /* Managing user's options. 
    options[0]: MinResSolve.
    options[1]: Plot Graph 1.
    options[2]: Plot Graph 2.
    */
    if (argc == 4){
        options[0] = atoi(argv[1]);
        options[1] = atoi(argv[2]);
        options[2] = atoi(argv[3]);
    } else {
        options[0] = 1;
        options[1] = 0;
        options[2] = 0;
    }


}


int main(int argc, char* argv[]){

    int options[3];
    input(options, argc, argv);

    if (options[0]){
        /* Testing the minimal residual linear solving function.
            [1 0 0]             [1]
        A = [0 1 0] + [1 1 1] * [1]
            [0 0 1]             [1]
        */
        cout << endl << "Testing ResMinSolve." << endl;
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
    }

    if (options[1]){
        cout << endl << "Plotting Graph 1." << endl;
        std::cout.precision(5);
        PlotGraph(10);
        cout << endl << "Graph 1 plotted." << endl << endl;
    }

    if (options[2]){
        cout << endl << "Plotting Graph 2." << endl;
        std::cout.precision(5);
        PlotGraph("matrix.txt");
        cout << endl << "Graph 2 plotted." << endl << endl;
    }

    return 0;
}
