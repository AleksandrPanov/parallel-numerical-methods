#pragma once
#include "CRS.h"

void generateMatrix(double *A, const int n, const int m)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A[i * m + j] = (std::rand() % 10) / 10.0 + 0.01;// -1.0;
}
void generateSimMatrix(double *A, const int n, const int m)
{
    for (int i = 0; i < n; i++)
        for (int j = i; j < m; j++)
        {
            A[i * m + j] = (std::rand() % 10) / 10.0 + 0.01;// -1.0;
            if (m > 1)
                A[j * m + i] = A[i * m + j];
        }
}
CRSMatrix generateCRSMatrix(const int n, const int len)
{
    CRSMatrix A = CRSMatrix();
    A.n = n;
    A.m = n;
    int count = 0;
    A.rowPtr.push_back(count);
    for (int i = 0; i < n; i++)
    {
        const int realLen = MIN(len, n - i);
        A.nz += realLen;
        count += realLen;
        if (realLen > 0)
        {
            A.rowPtr.push_back(count);
            for (int j = i; j < i + realLen; j++)
            {
                A.val.push_back((std::rand() % 10) / 10.0 + 0.01);
                A.colIndex.push_back(j);
            }
        }
    }
    return A;
}
void printMatrix(const vector<double> &A, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            cout << A[i*m + j] << " ";
        cout << "\n";
    }
    cout << "\n";
}
void printAsMatrix(const SLECRSMatrix &A) // нули, на месте где должны быть элементы симметричной матрицы
{
    int curIndx = 0;
    // идем по строкам:
    for (int i = 0; i < A.rowPtr.size() - 1; i++)
    {
        vector<double> row(A.m, 0.0);
        int rowElements = A.rowPtr[i + 1] - A.rowPtr[i];
        int endRow = curIndx + rowElements;

        for (curIndx; curIndx < endRow; curIndx++)
        {
            const int col = A.colIndexes[curIndx];
            row[col] = A.val[curIndx];
        }
        for (int j = 0; j < A.m; j++)
            cout << row[j] << " ";
        cout << "\n";
    }
    cout << "\n";
}
void printAsSimMatrix(const SLECRSMatrix &A) // вывод симметричной матрицы
{
    int curIndx = 0;
    vector<double> res(A.n*A.m, 0.0);
    // идем по строкам:
    for (int i = 0; i < A.n; i++)
    {
        const int rowElements = A.rowPtr[i + 1] - A.rowPtr[i];
        const int endRow = curIndx + rowElements;

        for (curIndx; curIndx < endRow; curIndx++)
        {
            const int col = A.colIndexes[curIndx];
            res[i*A.m + col] = A.val[curIndx];
            res[col*A.m + i] = A.val[curIndx];
        }
        for (int j = 0; j < A.m; j++)
            cout << res[i*A.m + j] << " ";
        cout << "\n";
    }
    cout << "\n";
}
double getError(const SLECRSMatrix &A, const vector<double> &vec, const vector<double> &x)
{
    const int n = x.size();
    vector<double> t(n * omp_get_max_threads());
    vector<double> res(n);
    //res = Ax
    A.mul(&x[0], &res[0], &t[0]);
    //res = Ax - b
    subtractVector(& res[0], &vec[0], 1.0, &res[0], n);
    return sqrt(scalar(&res[0], &res[0], n));
}