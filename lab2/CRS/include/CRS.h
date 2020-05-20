#include <omp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

#define MIN(a, b) ((a) <=(b) ? (a) : (b))

double scalar(const double *vec1, const double *vec2, int n)
{
    double res = 0.0;
    //#pragma omp parallel for reduction(+:res)
    for (int i = 0; i < n; i++)
    {
        res += vec1[i] * vec2[i];
    }
    return res;
}
void addVector(const double* vec1, const double* vec2, const double alpha, double* res, const int n)
{
    //#pragma omp parallel for
    //#pragma ivdep
    for (int i = 0; i < n; i++)
    {
        res[i] = vec1[i] + alpha * vec2[i];
    }
}
void subtractVector(const double* vec1, const double* vec2, const double alpha, double* res, const int n)
{
    //#pragma omp parallel for
    //#pragma ivdep
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
    CRSMatrix(vector<double> &A, int n): n(n), m(n)
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

    void printAsMatrix() const // нули, на месте где должны быть элементы симметричной матрицы
    {
        int curIndx = 0;
        // идем по строкам:
        for (int i = 0; i < rowPtr.size() - 1; i++)
        {
            vector<double> row(m, 0.0);
            int rowElements = rowPtr[i + 1] - rowPtr[i];
            int endRow = curIndx + rowElements;

            for (curIndx; curIndx < endRow; curIndx++)
            {
                const int col = colIndex[curIndx];
                row[col] = val[curIndx];
            }
            for (int j = 0; j < m; j++)
                cout << row[j] << " ";
            cout << "\n";
        }
        cout << "\n";
    }

    void printAsSimMatrix() const // вывод симметричной матрицы
    {
        int curIndx = 0;
        vector<double> res(n*m, 0.0);
        // идем по строкам:
        for (int i = 0; i < n; i++)
        {
            const int rowElements = rowPtr[i + 1] - rowPtr[i];
            const int endRow = curIndx + rowElements;

            for (curIndx; curIndx < endRow; curIndx++)
            {
                const int col = colIndex[curIndx];
                res[i*m + col] = val[curIndx];
                res[col*m + i] = val[curIndx];
            }
            for (int j = 0; j < m; j++)
                cout << res[i*m + j] << " ";
            cout << "\n";
        }
        cout << "\n";
    }

    void mul(const double* vec, double* res) const
    {
        //curIndx должен быть свой для каждого потока: indexThread*(n/numThreads)
        int curIndx = 0;
        //tmpTr должен быть свой для каждого потока
        vector<double> tmpTr(n);
        for (int i = 0; i < n; i++)
        {
            const int rowElements = rowPtr[i + 1] - rowPtr[i];
            const int endRow = curIndx + rowElements;
            res[i] = 0.0;
            for (curIndx; curIndx < endRow; curIndx++)
            {
                const int j = colIndex[curIndx];
                res[i] += val[curIndx] * vec[j];
                if (j != i)
                    tmpTr[j] += val[curIndx] * vec[i];
            }
            //должна быть редукция в res[i]
            res[i] += tmpTr[i];
            tmpTr[i] = 0.0;
        }
    }
};

void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
    /// initialization
    const int n = A.n;
    double alpha, beta;
    vector<double> r0(n), r1(n, 0.0), Ap(n, 0.0), p;

    //#pragma omp parallel for
    //#pragma ivdep
    for (int i = 0; i < n; i++)
        x[i] = 0.0;
    
    /// calculate
    // r0 = Ax0
    A.mul(x, &r0[0]);
    //#pragma omp parallel for
    //#pragma ivdep
    for (int i = 0; i < n; i++)
        r0[i] = b[i] - r0[i];

    //p = r0
    p = r0;

    for (count = 0; count < max_iter; count++)
    {
        //На softgader ошибка считается по другому
        double error = sqrt(scalar(&r0[0], &r0[0], n));
        if (error < eps)
            break;
        // alpha_i
        A.mul(&p[0], &Ap[0]);
        alpha = scalar(&r0[0], &r0[0], n) / scalar(&Ap[0], &p[0], n);
        
        // x_i+1
        addVector(&x[0], &p[0], alpha, &x[0], n);

        // r_i+1
        subtractVector(&r0[0], &Ap[0], alpha, &r1[0], n);

        // beta_i
        beta = scalar(&r1[0], &r1[0], n) / scalar(&r0[0], &r0[0], n);

        // p_i+1
        addVector(&r1[0], &p[0], beta, &p[0], n);

        // r_0 = r_i
        r0 = r1;
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