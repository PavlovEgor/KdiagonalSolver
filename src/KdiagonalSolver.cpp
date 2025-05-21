#include "KdiagonalSolver.hpp"
#include <cstdlib>
#include <iostream>


using std::vector;

KdiagonalSolver::KdiagonalSolver(
    const vector<vector<double>>& _diagonals, 
    const vector<int>& _indexes,
    const vector<double> _B)
{
    indexes     = std::move(_indexes);
    B           = _B;

    N = B.size();
    vector<double> zero_N(N, 0);
    
    diagonals.reserve(indexes.size());
    std::fill(diagonals.begin(), diagonals.end(), zero_N);

    k = 0;
    int ind_0;

    

    for (int i = 0; i < indexes.size(); i++)
    {   
    std::cout << "SEX" << std::endl;

        if (k < abs(indexes.at(i)))
        {
            k = abs(indexes.at(i));
        }
        
        if (indexes.at(i) < 0){
            for (int j=i; j < N; j++)
            {
                diagonals[i][j] += _diagonals[i][j-i];
            }
        } else if (indexes.at(i) == 0)
        {
            ind_0 = indexes.at(i);
        }
        
        else {
            for (int j=0; j < N-i; j++)
            {
                diagonals[i][j] += _diagonals[i][j];
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        diagonals[ind_0][i] = 1;
        B[i] /= _diagonals[ind_0][i];

        for (int j = 0; j < indexes.size(); j++)
        {
            diagonals[j][i] /= _diagonals[ind_0][i];
        }
    }

    
}

KdiagonalSolver::~KdiagonalSolver(){}

void KdiagonalSolver::solve(vector<double> x){
    int i, j, l, p;
    double frac, sum, tmp;

    vector<double> zero_k(k, 0);
    vector<vector<double>> P(N, zero_k);
    vector<double> R(N, 0);

    vector<vector<double>> P2(k+1, zero_k);
    vector<double> R2(k+1, 0);


    R[0] = B[0];
    for (i = 0; i < indexes.size(); i++)
    {
        if (indexes[i] > 0)
        {
            P[0][indexes[i]] = diagonals[i][0];
        }
    }

    for (i=1; i < N-1; i++)
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

        sum = 1;
        for (j = 0; j < indexes.size(); j++)
        {
            if (indexes[j] < 0)
            {
                sum += P2[-indexes[j]][0] * diagonals[j][i];
            }
        }

        frac /= sum;

        for (j = 0; j < indexes.size(); j++)
        {
            if (indexes[j] > 0){
                P[i][indexes[j]-1] = - diagonals[j][i] * frac;
            }
        }
        R[i] = B[i] * frac;
        
        for (j = 0; j < indexes.size(); j++)
        {
            if (indexes[j] > 0){
                for (l = 0; l < indexes.size(); l++)
                {
                    if (indexes[l] < 0)
                    {
                        P[i][indexes[j]-1] -= diagonals[l][i] * P2[-indexes[l]][indexes[j]] * frac;
                    }
                    
                }
                

            }
        }

        for ( j = 0; j < indexes.size(); j++)
        {
            if (indexes[j] < 0){
                R[i] -= diagonals[j][i] * R2[indexes[j]] * frac;
            }
        }
        

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

    sum = 1;
    for (j = 0; j < indexes.size(); j++)
    {
        if (indexes[j] < 0)
        {
            sum += P2[-indexes[j]][0] * diagonals[j][N-1];
        }
    }
    x[N-1] = B[N-1];

    for (j = 0; j < indexes.size(); j++)
    {
        if (indexes[j] < 0){
            x[N-1] -= diagonals[j][N-1] * R2[-indexes[j]];
        }
    }
    x[N-1] /= sum;

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
            x[i] += P[i][j] * x[i + 1 + j];
        }
    }
}