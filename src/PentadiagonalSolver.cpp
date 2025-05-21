#include <PentadiagonalSolver.hpp>


Pentadiagonal_matrix_solver::Pentadiagonal_matrix_solver(
    int _N,
    int _k,
    double *_first_up_diag,
    double *_second_up_diag,
    double *_first_down_diag,
    double *_second_down_diag)
{
    N = _N;
    k = _k;

    first_up_diag = _first_up_diag;
    second_up_diag = _second_up_diag;
    first_down_diag = _first_down_diag;
    second_down_diag = _second_down_diag;
}   

Pentadiagonal_matrix_solver::~Pentadiagonal_matrix_solver()
{

}

void Pentadiagonal_matrix_solver::pentadiagonal_solve(double *x, double *B){

    int i, j, l, p;
    double frac;
    double bi, ci, di, ei, tmp;

    double P[N][k];
    double R[N];

    double P2[k+1][k];
    double R2[k+1];

    R[0] = B[0];
    P[0][0] = -first_up_diag[0];
    P[0][k-1] = -second_up_diag[0];

    for (i = 1; i < k; i++)
    {
        di = first_down_diag[i];

        frac = 1 / (1 + P[i-1][0] * di);

        P[i][k-1] = - second_up_diag[i] * frac;
        for (j = 0; j < k-1; j++)
        {
            P[i][j] = - di * P[i-1][j+1] * frac;
        }
        P[i][0] += -first_up_diag[i] * frac;

        R[i] = (B[i] - di* R[i-1]) * frac;
        
    }

    for (i=k; i < N-1; i++)
    {
        
        for (j=0; j<k; j++)
        {
            P2[0][j] = 0;
            P2[1][j] = P[i-1][j];
        }

        P2[0][0] = 1;
        R2[0] = 0;
        R2[1] = R[i-1];

        for (j=2; j < k+1; j++)
        {
            for (l = 0; l < k - j + 1; l++)
            {
                P2[j][l] = P[i - j][j-1 + l];
            }
            for (l = k-j+1; l < k; l++)
            {
                P2[j][l] = 0;
            }
            
            R2[j] = R[i-j];

            for (l=1; l < j; l++)
            {
                tmp = P[i - j][j - l - 1];
                for (p=0; p < k; p++)
                {
                    P2[j][p] += tmp * P2[l][p];
                }
                R2[j] += tmp * R2[l];
            }
        }

        

        di = first_down_diag[i];
        ei = second_down_diag[i];

        frac = 1.0 / (1 + P[i-1][0] * di + P2[k][0] * ei);
        

        P[i][k-1] = - second_up_diag[i] * frac;

        
        for (j = 0; j < k-1; j++)
        {
            P[i][j] = - (di * P[i-1][j+1] + ei * P2[k][j+1]) * frac;
        }
        P[i][0] += -first_up_diag[i] * frac;

        R[i] = (B[i] - di* R[i-1] - ei * R2[k]) * frac;

    }

    for (j=0; j<k; j++)
    {
        P2[0][j] = 0;
        P2[1][j] = P[N-2][j];
    }

    P2[0][0] = 1;
    R2[0] = 0;
    R2[1] = R[N-2];

    for (j=2; j < k+1; j++)
    {
        for (l = 0; l < k - j + 1; l++)
        {
            P2[j][l] = P[N - 1 - j][j-1 + l];
        }
        for (l = k-j+1; l < k; l++)
        {
            P2[j][l] = 0;
        }
        
        R2[j] = R[N - 1 - j];

        for (l=1; l < j; l++)
        {
            tmp = P[N - 1 - j][j - l - 1];
            for (p=0; p < k; p++)
            {
                P2[j][p] += tmp * P2[l][p];
            }
            R2[j] += tmp * R2[l];
        }
    }


    di = first_down_diag[N-1];
    ei = second_down_diag[N-1];
    
    x[N-1] = (B[N-1] - di * R[N-2] - ei * R2[k]) / (1 + di * P[N-2][0] + ei * P2[k][0]);

    for (i = N-2; i > N - k - 1; i--)
    {
        x[i] = R[i];
        for (j=0; j<N - i - 1; j++)
        {
            x[i] += P[i][j] * x[i + 1 + j];
        }
    }

    for (i = N - k - 1; i >= 0; i--)
    {
        x[i] = R[i];
        for (j=0; j<k; j++)
        {
            x[i] += P[i][j] * x[i+ 1 + j];
        }
    }

}

