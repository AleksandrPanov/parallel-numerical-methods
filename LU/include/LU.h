#include <omp.h>
#include <vector>
#define indx(i, j, n) ((i)*(n)+(j))
#define MIN(a, b) ((a) <=(b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) < 0 ? -(a) : (a))

#include <iostream>
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

class BMatrix
{
    int firstIndex;
    double *A;
    int realSize;
    int n, m;
    int si, sj;
public:
    BMatrix(double *A, int rs, int n, int m, int si, int sj) : A(A), realSize(rs), n(n), m(m), si(si), sj(sj) 
    {
        firstIndex = indx(si, sj, realSize);
    };
    double get(int i, int j) const
    {
        return A[firstIndex + i * realSize + j];
    }
    double &operator() (int i, int j)
    {
        return A[firstIndex + i * realSize + j];
    }

    int row() const
    {
        return n;
    }
    int col() const
    {
        return m;
    }
    void print() const
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
                std::cout << (*this).get(i, j) << " ";
            std::cout << "\n";
        }
        std::cout << "\n";
    }

};

void blockMultMatrix(const BMatrix &A, const BMatrix &B, BMatrix &C);

void gaussU12(const BMatrix &L11, BMatrix &U12, const BMatrix &A12)
{
    //#pragma omp parallel for
    for (int r = 0; r < A12.row(); r++)
    {
        for (int i = U12.row() - 1; i >= 0; i--)
        {
            U12(i, r) = A12.get(i, r);
            #pragma ivdep
            for (int j = L11.col() - 1; j > i; j--)
                U12(i, r) -= L11.get(i, j) * U12(j, r);
        }
    }
}


void LU_Decomposition(double *A, double *L, double *U, int n)
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
        //этап 1, поиск L11 и U11
        for (int ved = 0; ved < n; ved++)
        {
            double el = U[indx(ved, ved, n)];
#pragma omp parallel for
            for (int i = ved + 1; i < n; i++)
            {
                double coeff = U[indx(i, ved, n)] / el;
                L[indx(i, ved, n)] = coeff;
                //#pragma ivdep
                for (int j = ved + 1; j < n; j++)
                    U[indx(i, j, n)] -= coeff * U[indx(ved, j, n)];

            }
        }
#pragma omp parallel for
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
            U[indx(i, j, n)] = 0.0;
}

void LU_Decomposition_block(double *A, double *L, double *U, int n) //block version
{
    const int bs = 2;
    for (int bi = 0; bi < n; bi += bs)
    {
        //этап 0, подготовка L11 и U11
        #pragma omp parallel for
        for (int i = bi; i < bi + bs; i++)
            L[indx(i, i, n)] = 1.0;

        #pragma omp parallel for
        for (int i = bi; i < bi + bs; i++)
            for (int j = i + 1; j < bi + bs; j++)
                L[indx(i, j, n)] = 0.0;

        #pragma omp parallel for
        for (int i = bi; i < bi + bs; i++)
            for (int j = bi; j < bi + bs; j++)
            U[indx(i, j, n)] = A[indx(i, j, n)];

        //этап 1, поиск L11 и U11
        for (int ved = bi; ved < MIN(bi + bs, n); ved++)
        {
            double el = U[indx(ved, ved, n)];
            #pragma omp parallel for
            for (int i = ved + 1; i < MIN(bi + bs, n); i++)
            {
                double coeff = U[indx(i, ved, n)] / el;
                L[indx(i, ved, n)] = coeff;
                #pragma ivdep
                for (int j = ved + 1; j < MIN(bi + bs, n); j++)
                    U[indx(i, j, n)] -= coeff * U[indx(ved, j, n)];

            }
        }     

        //этап 2 поиск U12
        BMatrix A11(A, n, bs, bs, bi, bi), L11(L, n, bs, bs, bi, bi), U11(U, n, bs, bs, bi, bi), 
            U12(U, n, MIN(n - bi, bs), MIN(n - bi - bs, bs), bi, bi + bs),
            A12(A, n, MIN(n - bi, bs), MIN(n - bi - bs, bs), bi, bi + bs);

        std::vector<double> res(bs*bs);
        A11.print();
        BMatrix R(&res[0], bs, bs, bs, 0, 0);
        blockMultMatrix(L11, U11, R);
        R.print();

        gaussU12(L11, U12, A12);

        //этап 3 поиск L21


        //этап 4 U -= L21 U12, не трогая U11
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

void blockMultMatrix(const BMatrix &A, const BMatrix &B, BMatrix &C)
{
    const int bs = 48;
#pragma omp parallel for
    for (int bi = 0; bi < A.row(); bi += bs)
        for (int bj = 0; bj < B.col(); bj += bs)
            for (int bk = 0; bk < A.col(); bk += bs)
                for (int i = bi; i < MIN(bi + bs, A.row()); i++)
                    for (int k = bk; k < MIN(bk + bs, B.col()); k++)
                        #pragma ivdep
                        for (int j = bj; j < MIN(bj + bs, A.col()); j++)
                        {
                            C(i, j) += A.get(i,k) * B.get(k,j);
                        }
}

double getMaxDiff(double *A, double *B, int n)
{
    double tmp = 0.0;
    for (int i = 0; i < n*n; i++)
    {
        tmp = MAX(tmp, ABS(A[i] - B[i]));
    }
    return tmp;
}