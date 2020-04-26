#include <omp.h>


#define indx(i, j, n) ((i)*(n)+(j))
#define MIN(a, b) ((a) <=(b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) < 0 ? -(a) : (a))


void LU_Decomposition(double *A, double *L, double *U, int n) //матрицы по строкам
{
    //  ||LU - A||/||A|| < 0.01
#pragma omp parallel for
    for (int i = 0; i < n; i++)
        L[indx(i, i, n)] = 1.0;

#pragma omp parallel for
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            L[indx(i, j, n)] = 0.0;


#pragma omp parallel for
    for (int i = 0; i < n*n; i++)
        U[i] = A[i];

    for (int ved = 0; ved < n; ved++)
    {
        double el = U[indx(ved, ved, n)];
        #pragma omp parallel for
        for (int i = ved + 1; i < n; i++)
        {
            double coeff = U[indx(i, ved, n)] / el;
            L[indx(i, ved, n)] = coeff;
            for (int j = ved + 1; j < n; j++)
                U[indx(i, j, n)] -= coeff * U[indx(ved, j, n)];

        }
    }
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
            U[indx(i, j, n)] = 0.0;
}

void multMatrix(double *A, double *B, double *C, int n)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
        for (int z = 0; z < n; z++)
#pragma ivdep
            for (int j = 0; j < n; j++)
            {
                C[indx(i, j, n)] += A[indx(i, z, n)] * B[indx(z, j, n)];
            }
}

void blockMultMatrix(double *A, double *B, double *C, int n)
{
    const int bs = 48;
    #pragma omp parallel for
    for (int bi = 0; bi < n; bi += bs)
        for (int bj = 0; bj < n; bj += bs)
            for (int bk = 0; bk < n; bk += bs)
                for (int i = bi; i < MIN(bi + bs, n); i++)
                    for (int k = bk; k < MIN(bk + bs, n); k++)
                        #pragma ivdep
                        for (int j = bj; j < MIN(bj + bs, n); j++)
                        {
                            C[i*n + j] += A[i*n + k] * B[k*n + j];
                        }
}

double getMaxDiff(double *A, double *B, int n)
{
    double tmp = 0.0;
    for (int i = 0; i < n*n; i++)
    {
        tmp = MAX(tmp, ABS(A[i] - B[i]));
    }
}