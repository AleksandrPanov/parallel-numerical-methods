#include <omp.h>
#include <vector>
using namespace std;

struct CRSMatrix
{
    int n; // ����� ����� � �������
    int m; // ����� �������� � �������
    int nz; // ����� ��������� ��������� � ����������� ������������ �������, ������� �� ���� ������� ���������
    vector<double> val; // ������ �������� ������� �� �������
    vector<int> colIndex; // ������ ������� ��������
    vector<int> rowPtr; // ������ �������� ������ �����
};

void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{

}