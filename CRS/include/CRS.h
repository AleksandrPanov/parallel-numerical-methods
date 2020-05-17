#include <omp.h>
#include <vector>
using namespace std;

struct CRSMatrix
{
    int n; // Число строк в матрице
    int m; // Число столбцов в матрице
    int nz; // Число ненулевых элементов в разреженной симметричной матрице, лежащих не ниже главной диагонали
    vector<double> val; // Массив значений матрицы по строкам
    vector<int> colIndex; // Массив номеров столбцов
    vector<int> rowPtr; // Массив индексов начала строк
};

void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{

}