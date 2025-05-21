#ifndef PENT_SOL
#define PENT_SOL

class Pentadiagonal_matrix_solver
{
private:
    int N; // size of matrix
    int k; // offset from the main diagonal


    double *first_up_diag;
    double *second_up_diag;
    double *first_down_diag;
    double *second_down_diag;


public:
    Pentadiagonal_matrix_solver(
        int _N, 
        int _k, 
        double *_first_up_diag, 
        double *_second_up_diag, 
        double *_first_down_diag, 
        double *_second_down_diag);
    ~Pentadiagonal_matrix_solver();

    void pentadiagonal_solve(double *x, double *B);
};

#endif