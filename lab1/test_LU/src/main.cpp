#include "LU.h"
void setZero(double *A, int n)
{
    for (int i = 0; i < n*n; i++)
        A[i] = 0.0;
}
int main()
{
    const int n = 20;
    double *A, *L, *U, *res;
    A = new double[n*n];
    L = new double[n*n];
    U = new double[n*n];
    res = new double[n*n];
    for (int i = 0; i < n*n; i++)
    {
        A[i] = i % (n / 3);
        res[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
    {
        A[indx(i, i, n)] = n;
    }
    std::cout << "A :\n";
    printMatrix(A, n);
    std::cout << "\n";

    //LU_Decomposition(A, L, U, n);
    //multMatrix(L, U, res, n);
    //printMatrix(L, n);
    //printMatrix(U, n);
    //printMatrix(res, n);
    //setZero(res, n);

    LU_Decomposition(A, L, U, n);
    //printMatrix(L, n);
    //printMatrix(U, n);
    blockMultMatrix(L, U, res, n, n, n, n);
    printMatrix(res, n);

    return 0;
}