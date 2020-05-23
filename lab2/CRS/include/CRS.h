#pragma once

#include <omp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

#define MIN(a, b) ((a) <=(b) ? (a) : (b))

double scalar(const double *vec1, const double *vec2, const int n)
{
    double res = 0.0;
    #pragma omp parallel for reduction(+:res)
    for (int i = 0; i < n; i++)
    {
        res += vec1[i] * vec2[i];
    }
    return res;
}
double gerError(const double* x0, const double* x1, const double* b, const int n)
{
    double error1 = 0.0;
    double error2 = 0.0;
    #pragma omp parallel for reduction(+:error1, error2)
    for (int i = 0; i < n; i++)
    {
        error1 += (x1[i] - x0[i]) * (x1[i] - x0[i]);
        error2 += b[i] * b[i];
    }
    return sqrt(error1) / sqrt(error2);
}
void copy_ar(const double *vec, double *res, int n)
{
    #pragma omp parallel for
    #pragma ivdep
    for (int i = 0; i < n; i++)
    {
        res[i] = vec[i];
    }
}
void addVector(const double* vec1, const double* vec2, const double alpha, double* res, const int n)
{
    #pragma omp parallel for
    #pragma ivdep
    for (int i = 0; i < n; i++)
    {
        res[i] = vec1[i] + alpha * vec2[i];
    }
}
void subtractVector(const double* vec1, const double* vec2, const double alpha, double* res, const int n)
{
    #pragma omp parallel for
    #pragma ivdep
    for (int i = 0; i < n; i++)
    {
        res[i] = vec1[i] - alpha * vec2[i];
    }
}

struct CRSMatrix
{
    int n; // Число строк в матрице
    int m; // Число столбцов в матрице
    int nz; // Число ненулевых элементов в разреженной симметричной матрице, лежащих не ниже главной диагонали
    vector<double> val; // Массив значений матрицы по строкам
    vector<int> colIndex; // Массив номеров столбцов
    vector<int> rowPtr; // Массив индексов начала строк
    CRSMatrix():n(0), m(0), nz(0){}
    CRSMatrix(vector<double> &A, int n) : n(n), m(n)
    {
        rowPtr.push_back(0);
        for (int i = 0; i < n; i++)
        {
            int endRow = rowPtr[i];
            // считаем, что на вход идут симметричные матрицы
            for (int j = i; j < m; j++)
            {
                double tmp = A[i*m + j];
                if (tmp != 0.0)
                {
                    val.push_back(tmp);
                    colIndex.push_back(j);
                    endRow++;
                }
            }
            rowPtr.push_back(endRow);
        }
        nz = rowPtr.back();
    }
};

struct SLECRSMatrix
{
    const int &n; // Число строк в матрице
    const int &m; // Число столбцов в матрице
    const int &nz; // Число ненулевых элементов в разреженной симметричной матрице, лежащих не ниже главной диагонали
    const vector<double> &val; // Массив значений матрицы по строкам
    const vector<int> &colIndexes; // Массив номеров столбцов
    const vector<int> &rowPtr; // Массив индексов начала строк
    SLECRSMatrix(const CRSMatrix &matr):n(matr.n), m(matr.m), nz(matr.nz), val(matr.val), colIndexes(matr.colIndex), rowPtr(matr.rowPtr){}
    void mul(const double* vec, double* res, double* t) const
    {
        const int numThreads = omp_get_max_threads();
        #pragma omp parallel
        {
            const int indxThread = omp_get_thread_num();
            const int size = n / numThreads;
            const int start = indxThread * size;
            int ostatok = 0;
            if (indxThread == numThreads - 1)
                ostatok = n % numThreads;
            const int end = start + size + ostatok;
            double *tmp = t + indxThread * n;
            int elIndx = rowPtr[start];
            for (int i = start; i < end; i++)
            {
                const int endRow = elIndx + (rowPtr[i + 1] - rowPtr[i]);
                res[i] = 0.0;
                for (elIndx; elIndx < endRow; elIndx++)
                {
                    const int j = colIndexes[elIndx];
                    const double vv = val[elIndx];
                    res[i] += vv * vec[j];
                    if (j != i)
                        tmp[j] += vv * vec[i];
                }
            }
        }
        //редукция в res[i]
        for (int thr = numThreads - 1; thr >= 0; thr--)
        {
            double *tmp = t + thr * n;
            #pragma omp parallel for
            #pragma ivdep
            for (int i = 0; i < n; i++)
            {
                res[i] += tmp[i];
                tmp[i] = 0.0;
            }
        }
    }
};
void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
    /// initialization
    const int n = A.n;
    double alpha, beta;
    SLECRSMatrix mA(A);
    vector<double> r0(n, 0.0), r1(n, 0.0), Ap(n, 0.0), p(n, 0.0), x0(n, 0.0), tmp(omp_get_max_threads() * n, 0.0);
    #pragma omp parallel for
    #pragma ivdep
    for (int i = 0; i < n; i++)
        x[i] = 0.0;
    
    /// calculate
    // r0 = Ax0
    mA.mul(x, &r0[0], &tmp[0]);
    // r0 = b - Ax0
    #pragma omp parallel for
    #pragma ivdep
    for (int i = 0; i < n; i++)
        r0[i] = b[i] - r0[i];

    //p = r0
    copy_ar(&r0[0], &p[0], n);

    for (count = 0; count < max_iter; count++)
    {
        //На softgader ошибка считается по другому
        //double error = sqrt(scalar(&r0[0], &r0[0], n));
        //if (error < eps)
        //    break;

        // alpha_i
        mA.mul(&p[0], &Ap[0], &tmp[0]);
        alpha = scalar(&r0[0], &r0[0], n) / scalar(&Ap[0], &p[0], n);
        
        // копируем в x0 "предыдущий" ответ
        // x_i = x_i+1
        copy_ar(&x[0], &x0[0], n);

        // x_i+1
        addVector(&x[0], &p[0], alpha, &x[0], n);

        // r_i+1
        subtractVector(&r0[0], &Ap[0], alpha, &r1[0], n);

        // beta_i
        beta = scalar(&r1[0], &r1[0], n) / scalar(&r0[0], &r0[0], n);

        // p_i+1
        addVector(&r1[0], &p[0], beta, &p[0], n);

        // r_0 = r_i
        copy_ar(&r1[0], &r0[0], n);

        double error = gerError(&x0[0], x, b, n);
        if (error < eps)
        {
            count++;
            return;
        }
    }

}

void blockMultMatrix(double *A, double *B, double *C, int n1, int m1, int n2, int m2)
{
    const int bs = 48;
#pragma omp parallel for
    for (int bi = 0; bi < n1; bi += bs)
        for (int bj = 0; bj < m2; bj += bs)
            for (int bk = 0; bk < m1; bk += bs)
                for (int i = bi; i < MIN(bi + bs, n1); i++)
                    for (int k = bk; k < MIN(bk + bs, m1); k++)
                    {
                        const int bjmin = MIN(bj + bs, m2);
#pragma ivdep
                        for (int j = bj; j < bjmin; j++)
                        {
                            C[i*m2 + j] += A[i*m1 + k] * B[k*m2 + j];
                        }
                    }
}