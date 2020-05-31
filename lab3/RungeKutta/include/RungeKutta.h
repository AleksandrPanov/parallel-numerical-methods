#pragma once

#include <omp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

#define MIN(a, b) ((a) <=(b) ? (a) : (b))

class heat_task {
public:
    double T; // ������ �������, � ������� ���������� ���������������� u(x, t)
    double L; // ����� �������
    int n; // ������ ����� �� x
    int m; // ������ ����� �� t
    double initial_condition(double x); // �������, �������� ��������� �������
    double left_condition(double t); // �������, �������� ��������� ������� ��� x = 0
    double right_condition(double t); // �������, �������� ��������� ������� ��� x = L
    double f(double x, double t); // �������, �������� ������� �����������
};
void heat_equation_runge_kutta(heat_task task, double * v)
{

}
// ������� �������� � ���������� ��������� ���������� :
// task � ���������� ������ ���� heat_task, � ������� �������� �������� ������.
// ������������� ���������� ������� �� ���������� ������ �� ������� ������.
