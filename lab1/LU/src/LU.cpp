#include "LU.h"
#include <cstdlib>
#include <ctime>
#include <vector>
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
    //omp_set_num_threads(8);
    const int n = 4000;
    double *A, *L, *U, *res;
    A = new double[n*n];
    L = new double[n*n];
    U = new double[n*n];
    res = new double[n*n];
    #pragma omp parallel for
    for (int i = 0; i < n*n; i++) 
    {
        A[i] = 0.0;
        res[i] = 0.0;
    }

    generateMatrix(L, n);
    //printMatrix(L, n);
    transposeMatrix(L, U, n);

    blockMultMatrix(L, U, A, n, n, n, n);

    std::vector<double> acopy(A, A + n*n);

   //printMatrix(A, n);

    auto t1 = clock();
    LU_Decomposition(A, L, U, n);
    auto t2 = clock() - t1;
    std::cout << "time " << t2 << "\n";
    //printMatrix(L, n);
    //printMatrix(U, n);


    blockMultMatrix(L, U, res, n, n, n, n);
    //printMatrix(res, n);
    std::cout << getMaxDiff(res, &acopy[0], n);
    return 0;
}