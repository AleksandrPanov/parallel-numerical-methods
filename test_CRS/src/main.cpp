#include "CRS.h"
#include <ctime>

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
    vector<double> Atmp = { 1, 2, 3, 2, 1, 4, 3, 4, 1 }, vec = {1, -1, 2};
    CRSMatrix A = CRSMatrix(Atmp, 3);
    A.printAsSimMatrix();

    print(A.mul(vec), 3, 1);
    return 0;
}