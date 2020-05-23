#include "CRS.h"
#include "print_generate.h"
#include <cstdlib>
#include <ctime>


int main()
{
    // 1 thread 20
    // 2 thread 10.8
    // 4 thread 6.2
    // 6 thread 5.1
    // 8 thread 4.7
    omp_set_num_threads(4);
    const int n = 18000;
    const int len = 3;
    CRSMatrix tmp = generateCRSMatrix(n, len);
    SLECRSMatrix A(tmp);
    vector<double> vec(n), res(n);
    generateMatrix(&vec[0], n, 1);
    
    //printMatrix(vec, n, 1);
    //printAsSimMatrix(A);

    int count = 0;
    int time1 = clock();
    SLE_Solver_CRS(tmp, &vec[0], 0.0, 80*1000, &res[0], count);
    int time2 = clock() - time1;
    cout << "time: " << time2 << "\n";

    //printMatrix(res, n, 1);

    std::cout << "num iteration: " << count << "\n";
    std::cout << "error: " << getError(A, vec, res) << "\n";
    return 0;
}