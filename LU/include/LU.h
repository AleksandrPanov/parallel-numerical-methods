#include <omp.h>
#define indx(i, j, n) ((i)*(n)+(j))
#define MIN(a, b) ((a) <=(b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) < 0 ? -(a) : (a))

class BMatrix
{
    const int firstIndex;
    double *A;
    const int realSize;
    const int n, m;
    const int si, sj;
public:
    BMatrix(double *A, int rs, int n, int m, int si, int sj) : A(A), realSize(rs), n(n), m(m), si(si), sj(sj), firstIndex(indx(si, sj, rs)) {};
    __forceinline double &operator() (int i, int j)
    {
        return A[firstIndex + i * realSize + j];
    }
    __forceinline const double &operator() (int i, int j) const
    {
        return A[firstIndex + i * realSize + j];
    }
    __forceinline int row() const
    {
        return n;
    }
    __forceinline int col() const
    {
        return m;
    }
    int end_row() const
    {
        return row() + si;
    }
    int end_col() const
    {
        return col() + sj;
    }
    //void print() const
    //{
    //    for (int i = 0; i < n; i++)
    //    {
    //        for (int j = 0; j < m; j++)
    //            std::cout << (*this)(i, j) << " ";
    //        std::cout << "\n";
    //    }
    //    std::cout << "\n";
    //}
    void blockMultMatrix(const BMatrix &mA, const BMatrix &mB, const int bs)
    {
        #pragma omp parallel for
        for (int bi = 0; bi < mA.row(); bi += bs)
            for (int bj = 0; bj < mB.col(); bj += bs)
                for (int bk = 0; bk < mA.col(); bk += bs)
                    for (int i = bi; i < MIN(bi + bs, mA.row()); i++)
                    {
                        const double *A = mA.A + mA.firstIndex + i * realSize;
                        double *mThis = this->A + firstIndex + i * realSize;
                        const int bkmin = MIN(bk + bs, mA.col());
                        for (int k = bk; k < bkmin; k++)
                        {
                            const double *B = mB.A + mB.firstIndex + k * realSize;
                            const int bjmin = MIN(bj + bs, mB.col());
                            #pragma ivdep
                            for (int j = bj; j < bjmin; j++)
                                mThis[j] -= A[k] * B[j]; //A[indxTmp + j] -= mA(i, k) * B(k, j);
                        }
                    }
    }
};

void LU(BMatrix &L11, BMatrix &U11, const BMatrix &A11)
{
    #pragma omp parallel for
    for (int i = 0; i < L11.row(); i++)
        for (int j = 0; j < L11.col(); j++)
        {
            if (i != j)
                L11(i, j) = 0.0;
            else L11(i, j) = 1.0;
        }

    #pragma omp parallel for
    for (int i = 0; i < A11.row(); i++)
        for (int j = 0; j < A11.col(); j++)
            U11(i,j) = A11(i,j);

    for (int ved = 0; ved < U11.row(); ved++)
    {
        double el = U11(ved, ved);
        #pragma omp parallel for
        for (int i = ved + 1; i < U11.row(); i++)
        {
            double coeff = U11(i, ved) / el;
            L11(i, ved) = coeff;
            #pragma ivdep
            for (int j = 0; j <= ved; j++)
                U11(i, j) = 0.0;
            #pragma ivdep
            for (int j = ved + 1; j < U11.col(); j++)
                U11(i, j) -= coeff * U11(ved, j);
        }
    }
}
void gaussU12(const BMatrix &L11, BMatrix &U12, const BMatrix &A12)
{
    #pragma omp parallel for
    for (int r = 0; r < A12.col(); r++)
    {
        #pragma ivdep
        for (int i = 0; i < U12.row(); i++)
            U12(i, r) = A12(i, r);

        for (int i = 1; i < U12.row(); i++)
        {
            #pragma novector
            for (int j = 0; j < i; j++)
                U12(i, r) -= L11(i, j) * U12(j, r);
        }
    }
}
void gaussL21(BMatrix &L21, const BMatrix& U11, const BMatrix& A21)
{
#pragma omp parallel for
    for (int rr = 0; rr < A21.row(); rr++)
    {
        #pragma ivdep
        for (int j = 0; j < L21.col(); j++)
            L21(rr, j) = A21(rr, j);

        for (int j = 0; j < L21.col(); j++)
        {
            #pragma novector
            for (int k = 0; k < j; k++)
                L21(rr, j) -= L21(rr, k) * U11(k, j);
            L21(rr, j) /= U11(j, j);
        }
    }
}
void LU_Decomposition(double *A, double *L, double *U, int n) //block version
{
    const int bs = 32;
    for (int bi = 0; bi < n; bi += bs)
    {
        //этап 0, подготовка L11 и U11
        BMatrix A11(A, n, MIN(n - bi, bs), MIN(n - bi, bs), bi, bi),
                L11(L, n, MIN(n - bi, bs), MIN(n - bi, bs), bi, bi),
                U11(U, n, MIN(n - bi, bs), MIN(n - bi, bs), bi, bi);
        //этап 1, поиск L11 и U11
        LU(L11, U11, A11);

        //этап 2 поиск U12
        BMatrix U12(U, n, MIN(n - bi, bs), n - U11.end_col(), bi, bi + bs),
                A12(A, n, MIN(n - bi, bs), n - A11.end_col(), bi, bi + bs);

        gaussU12(L11, U12, A12);

        //этап 3 поиск L21
        BMatrix L21(L, n, n - L11.end_row(), MIN(n - bi, bs), bi + bs, bi),
                A21(A, n, n - A11.end_row(), MIN(n - bi, bs), bi + bs, bi);
        gaussL21(L21, U11, A21);

        //этап 4 A -= L21 U12
        BMatrix A22(A, n, n - A11.end_row(), n - A11.end_col(), bi + bs, bi + bs);
        A22.blockMultMatrix(L21, U12, 32);
    }
    BMatrix Ures(U, n, n, n, 0, 0), Lres(L, n, n, n, 0, 0);
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
            Ures(i, j) = 0.0;
        for (int j = i + 1; j < n; j++)
            Lres(i, j) = 0.0;
    }
}

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
    return tmp;
}
void LU_Decomposition_notBlock(double *A, double *L, double *U, int n)
{
    // ||LU - A||/||A|| < 0.01
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
            #pragma ivdep
            for (int j = ved + 1; j < n; j++)
                U[indx(i, j, n)] -= coeff * U[indx(ved, j, n)];

        }
    }
#pragma omp parallel for
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
            U[indx(i, j, n)] = 0.0;
}