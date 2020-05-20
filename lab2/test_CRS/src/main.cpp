 #include "CRS.h"
#include <ctime>

void generateMatrix(double *A, const int n, const int m)
{
    for (int i = 0; i < n; i++)
        for (int j = i; j < m; j++)
        {
            A[i * m + j] = (std::rand() % 21) / 10.0 - 1.0;
            if (m > 1)
                A[j * m + i] = A[i * m + j];
        }
}

void print(const vector<double> &A, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            cout << A[i*m + j] << " ";
        cout << "\n";
    }
    cout << "\n";
}
int main()
{
    const int n = 3;
    vector<double> Atmp = vector<double>(n*n);
    Atmp = { 1, 2, 3, 2, 1, 4, 3, 4, 1 };
    vector<double> vec = vector<double>(n);
    vec = {1, -1, 2};
    //generateMatrix(&Atmp[0], n, n);
    //generateMatrix(&vec[0], n, 1);

    CRSMatrix A = CRSMatrix(Atmp, n);
    A.printAsSimMatrix();

    vector<double> res(n);
    vector<double> res1(n);
    vector<double> res2(n);

    
    blockMultMatrix(&Atmp[0], &vec[0], &res1[0], n, n, n, 1);
    print(res1, n, 1);

    A.mul(&vec[0], &res2[0]);
    print(res2, n, 1);

    int count = 0;
    SLE_Solver_CRS(A, &vec[0], 0.001, 100, &res[0], count);
    print(res, n, 1);

    std::cout << "num iteration: " << count;
    return 0;
}