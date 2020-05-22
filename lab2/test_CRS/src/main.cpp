#include "CRS.h"
#include "print_generate.h"
#include <ctime>

int main()
{
    const int n = 10;
    vector<double> Atmp = vector<double>(n*n);
    //Atmp = { 1, 2, 3, 2, 1, 4, 3, 4, 1 };
    vector<double> vec = vector<double>(n);
    //vec = {1, -1, 2};
    generateSimMatrix(&Atmp[0], n, n);
    generateMatrix(&vec[0], n, 1);
    vector<double> omp_tmp(omp_get_max_threads() * n, 0.0);


    CRSMatrix tmp(CRSMatrix(Atmp, n));
    SLECRSMatrix A(tmp);
    printAsSimMatrix(A);

    vector<double> res(n);
    vector<double> res1(n);
    vector<double> res2(n);

    
    //blockMultMatrix(&Atmp[0], &vec[0], &res1[0], n, n, n, 1);
    //printMatrix(res1, n, 1);
    //
    //A.mul(&vec[0], &res2[0], &omp_tmp[0]);
    //printMatrix(res2, n, 1);

    int count = 0;
    int time1 = clock();
    SLE_Solver_CRS(tmp, &vec[0], 0.01, 10000, &res[0], count);
    int time2 = clock() - time1;
    cout << "time: " << time2 << "\n";
    printMatrix(res, n, 1);

    std::cout << "num iteration: " << count;
    return 0;
}