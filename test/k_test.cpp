#include <iostream>
#include <KdiagonalSolver.hpp>

using std::vector;

int main(){
    int N = 10;
    int k = 3;

    vector<double> sud = { 0.56838519, 0.59106732, 
                    0.60376811, 0.06784427, 
                    0.97604943, 0.68617165, 
                    0.03754737};

    vector<double> fud = { 0.21496537, 0.23166263, 
                    0.55616111, 0.22583128, 
                    0.10656655, 0.54675014,
                    0.37189687, 0.38990937, 
                    0.98233873};

    vector<double> fdd = {0.68343727, 
                    0.89475701, 0.5148094, 
                    0.03900925, 0.68854438,
                    0.89352028, 0.2643208,
                    0.26646279, 0.9347437 };

    vector<double> sdd = {0.88109126,
                    0.62575212, 0.13503155,
                    0.30101654, 0.21706122, 
                    0.4206771,  0.10920916};

    vector<double> md = {1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1};

    const vector<vector<double>> _matrix = {md, sud, fud, fdd, sdd};

    const vector<int> ind = {0, 3, 1, -1, -3};


    vector<double> x(N);
    vector<double> B = {0.94122422, 0.86527486, 0.66331387, 0.51497868,
                    0.69172562, 0.99599477, 0.05788199, 0.26399862,
                    0.39642472, 0.60151007};
                    
    KdiagonalSolver test_LS = KdiagonalSolver(_matrix, ind, B);

    test_LS.solve(x);

    std::cout.precision(9);

    for ( int i = 0; i < 10; i++)
    {
        std::cout << x[i] <<" ";
    }
    std::cout << std::endl; 

    std::cout << "-1.06064525 2.14129119 -2.47122364 2.71218541 0.03613519 -0.48000190 0.02015260 -0.75710842 2.58505392 -1.81705365" << std::endl;

    return 0;
}