#include "KdiagonalSolver.hpp"
#include <cstdlib>
#include <iostream>
#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end

using std::vector;

KdiagonalSolver::KdiagonalSolver(
    const vector<vector<double>>& diagonals, 
    const vector<int>& indexes,
    const vector<double> _B) :
    B(_B)
{
    int i, j, count_of_diag;
    vector<double> Main_diagonal(N);

    count_of_diag = indexes.size();
    N = B.size();
    k_low = -*std::min_element(std::begin(indexes), std::end(indexes));
    k_up = *std::max_element(std::begin(indexes), std::end(indexes));

    Lower_diagonals.resize(k_low, vector<double>(N, 0.0));
    Upper_diagonals.resize(k_up, vector<double>(N, 0.0));

    for (j = 0; j < count_of_diag; j++)
    {
        if (indexes[j] < 0){

            for (i=-indexes[j]; i < N; i++)
            {
                Lower_diagonals[-indexes[j] - 1][i] += diagonals[j][i + indexes[j]];
            }

        } else if (indexes[j] == 0)
        {
            Main_diagonal = diagonals[j];
        }
        
        else {
            for (i=0; i < N-indexes[j]; i++)
            {
                Upper_diagonals[indexes[j] - 1][i] += diagonals[j][i];
            }
        }
    }
    

    for (i = 0; i < N; i++)
    {
        B[i] /= Main_diagonal[i];

        for (j = 0; j < k_low; j++)
        {
            Lower_diagonals[j][i] /= Main_diagonal[i];
        }
        for (j = 0; j < k_up; j++)
        {
            Upper_diagonals[j][i] /= Main_diagonal[i];
        }
    }    
}

KdiagonalSolver::~KdiagonalSolver(){
}

void KdiagonalSolver::solve(vector<double>& x){
    int i, j, l, p;
    double frac, sum, tmp;

    vector<vector<double>> P(N, vector<double>(k_up, 0.0));
    vector<double> R(N, 0);

    vector<vector<double>> Q(k_low+1, vector<double>(k_up, 0.0));
    vector<double> W(k_low+1, 0);

    R[0] = B[0];
    for (l = 0; l < k_up; l++){
        P[0][l] = -Upper_diagonals[l][0];
    }


    for (i=1; i < N; i++)
    {
        Q[0][0] = 1;
        W[0] = 0;
        
        for (j=0; j<k_up; j++)
        {
            Q[0][j] = 0;
            Q[1][j] = P[i-1][j];
        }
        W[1] = R[i-1];

        for (l=2; l <= k_low; l++)
        {
            if (i - l < 0){
                break;
            }

            W[l] = R[i - l];

            for (j = 0; j < k_up; j++)
            {
                Q[l][j] = 0;
            }

            Q[l][0] = P[i-l][l - 1];

            for (j = 0; j < k_up - l; j++)
            {
                Q[l][j + 1] = P[i-l][l + j];
            }

            for (j = 0; j <= l - 2; j++)
            {
                W[l] += P[i-l][j] * W[l-j-1];

                for (p = 0; p < k_up; p++)
                {
                    Q[l][p] += Q[l - j - 1][p] * P[i - l][j];
                }
            }            
        }

        sum = 1;
        for (l = 0; l < k_low; l++)
        {
            sum += Lower_diagonals[l][i] * Q[l+1][0];
        }
        frac = 1 / sum;

        if (i != N-1) {
        sum = B[i];
        for (l = 0; l < k_up; l++)
        {
            sum -= Lower_diagonals[l][i] * W[l + 1];
        }

        R[i] = sum * frac;

        for (j = 0; j < k_up; j++)
        {
            P[i][j] -= Upper_diagonals[j][i];
        }

        for (j = 0; j < k_up-1; j++)
        {
            for (l = 0; l < k_low; l++)
            {
                P[i][j] -= Lower_diagonals[l][i] * Q[l+1][j+1];
            }

            P[i][j] *= frac;
        }
        P[i][k_up-1] *= frac;
        
    } else if (i == N-1)
    {
        x[N-1] = B[N-1];
        for (l = 0; l < k_low; l++)
        {
            x[N-1] -= Lower_diagonals[l][N-1] * W[l+1];
        }
        x[N-1] *= frac;
    }
    }


    for (i = N-2; i > N - k_up - 1; i--)
    {
        x[i] = R[i];
        for (j=0; j < N - i - 1; j++)
        {
            x[i] += P[i][j] * x[i + 1 + j];
        }
    }

    for (i = N - k_up - 1; i >= 0; i--)
    {
        x[i] = R[i];
        for (j=0; j<k_up; j++)
        {
            x[i] += P[i][j] * x[i + 1 + j];
        }
    }
}