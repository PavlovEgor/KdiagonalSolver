#ifndef K_SOL
#define K_SOL

#include <memory>
#include <vector>

class KdiagonalSolver
{
private:
    int N; // size of matrix
    int k_low; // offset from the main diagonal
    int k_up; // offset from the main diagonal

    std::vector<std::vector<double>> Upper_diagonals;
    std::vector<std::vector<double>> Lower_diagonals;

    std::vector<double> B;


public:
    KdiagonalSolver(
        const std::vector<std::vector<double>> &_diagonals,
        const std::vector<int> &_indexes,
        const std::vector<double> _B);

    ~KdiagonalSolver();

    void solve(std::vector<double>& x);
};

#endif