#include "LU.h"
#include <cstdlib>
#include <ctime>


int main()
{
    const int n = 4;
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
    printMatrix(A, n);
    LU_Decomposition_block(A, L, U, n);
    //auto t3 = clock() - t2;
    //std::cout << t3 << "\n";
    //
    //multMatrix(L, U, res, n);
    ////printMatrix(res, n);
    //std::cout << getMaxDiff(res, A, n);
    return 0;
}