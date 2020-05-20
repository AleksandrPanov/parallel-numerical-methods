#include <omp.h>
#include <vector>
#include <iostream>
using namespace std;

struct CRSMatrix
{
    int n; // ����� ����� � �������
    int m; // ����� �������� � �������
    int nz; // ����� ��������� ��������� � ����������� ������������ �������, ������� �� ���� ������� ���������
    vector<double> val; // ������ �������� ������� �� �������
    vector<int> colIndex; // ������ ������� ��������
    vector<int> rowPtr; // ������ �������� ������ �����
    CRSMatrix(vector<double> &A, int n): n(n), m(n)
    {
        rowPtr.push_back(0);
        for (int i = 0; i < n; i++)
        {
            int endRow = rowPtr[i];
            // �������, ��� �� ���� ���� ������������ �������
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

    void printAsMatrix() const // ����, �� ����� ��� ������ ���� �������� ������������ �������
    {
        int curIndx = 0;
        // ���� �� �������:
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

    void printAsSimMatrix() const // ����� ������������ �������
    {
        int curIndx = 0;
        vector<double> res(n*m, 0.0);
        // ���� �� �������:
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

    vector<double> mul(const vector<double> &vec)
    {
        vector<double> res(n);
        //curIndx ������ ���� ���� ��� ������� ������: indexThread*(n/numThreads)
        int curIndx = 0;
        vector<double> tmpTr(n);
        for (int i = 0; i < n; i++)
        {
            const int rowElements = rowPtr[i + 1] - rowPtr[i];
            const int endRow = curIndx + rowElements;
            for (curIndx; curIndx < endRow; curIndx++)
            {
                const int j = colIndex[curIndx];
                res[i] += val[curIndx] * vec[j];
                if (j != i)
                    tmpTr[j] += val[curIndx] * vec[i];
            }
            res[i] += tmpTr[i];
        }
        
        return res;
    }
};

void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{

}