#include "LU.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

void printMatrix(double *A, int n)
{
    for (int i = 0; i < n*n; i++)
    {
        if (i % n == 0 && i != 0)
            std::cout << "\n";
        std::cout << A[i] << " ";
    }
    std::cout << "\n";
    std::cout << "\n";
}
void generateMatrix(double *A, int n)
{
    for (int i = 0; i < n*n; i++)
    {
        A[i] = (std::rand() % 21) / 10.0 - 1.0;
    }
}
void transposeMatrix(double *A, double *B, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            B[indx(j, i, n)] = A[indx(i, j, n)];
}
int main()
{
    const int n = 2000;
    double *A, *L, *U, *res;
    A = new double[n*n];
    L = new double[n*n];
    U = new double[n*n];
    res = new double[n*n];
    for (int i = 0; i < n*n; i++) 
    {
        A[i] = 0.0;
        res[i] = 0.0;
    }

    generateMatrix(L, n);
    transposeMatrix(L, U, n);

    auto t1 = clock();
    blockMultMatrix(L, U, A, n);
    auto t2 = clock() - t1;
    std::cout << t2 << "\n";
    //multMatrix(L, U, A, n);
    //auto t3 = clock() - t2;
    //std::cout << t3 << "\n";
    //std::cout << t2/(double)t3 << "\n";
    //printMatrix(A, n);
    LU_Decomposition(A, L, U, n);
    auto t3 = clock() - t2;
    std::cout << t3 << "\n";

    multMatrix(L, U, res, n);
    //printMatrix(res, n);
    std::cout << getMaxDiff(res, A, n);
    //printMatrix(U, n);
    //LU_Decomposition(A, L, U, n);

    return 0;
}