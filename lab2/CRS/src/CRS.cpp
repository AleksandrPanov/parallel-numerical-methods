#include "CRS.h"
#include "print_generate.h"
#include <cstdlib>
#include <ctime>


int main()
{
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
    SLE_Solver_CRS(tmp, &vec[0], 0.001 * 0.001, 1000*1000, &res[0], count);
    int time2 = clock() - time1;
    cout << "time: " << time2 << "\n";

    //printMatrix(res, n, 1);

    std::cout << "num iteration: " << count << "\n";
    std::cout << "error: " << getError(A, vec, res) << "\n";
    return 0;
}