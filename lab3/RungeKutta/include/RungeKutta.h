#pragma once
#include <algorithm>
class heat_task {
public:
    int n = 50000; // размер сетки по x
    int m = 50000; // размер сетки по t
    double L = 100; // длина стержн€
    double T = 0.1; // момент времени, в который необходимо аппроксимировать u(x, t)
    double initial_condition(double x) { // функци€, задающа€ начальное условие
        return pow(x, 4);
    }
    double left_condition(double t) { // функци€, задающа€ граничное условие при x = 0
        return pow(t, 4);
    }
    double right_condition(double t) { // функци€, задающа€ граничное условие при x = L
        return pow(L, 4) + pow(t, 4);
    }
    double f(double x, double t) { // функци€, задающа€ внешнее воздействие
        return 4 * pow(t, 3) - 12 * pow(x, 2);
    }
    double u(double x, double t) {
        return pow(x, 4) + pow(t, 4);
    }
    double error(double v, double x, double t) {
        return abs(v - u(x, t));
    }
};

#include <algorithm>
#include <omp.h>
#include <vector>
#include <iostream>
using namespace std;

class heat_task;

void calc_k1(heat_task &task, const double *v, double *k1, double t)
{
    const int n = task.n;
    const double h = task.L / n;
    const double coeff = 1 / (h*h);
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        k1[i] = coeff*(v[i + 1] - 2.0 * v[i] + v[i - 1]) + task.f(h*i, t);
    }
}
void calc_k2(heat_task &task, const double *v, const double *k1, double *k2, double t, double *tmp)
{
    const int n = task.n;
    const double h = task.L / n;
    const double coeff = 1 / (h*h);
    const double dt = task.T/(task.m * 2.0);
    tmp[0] = task.left_condition(t);
    tmp[n] = task.right_condition(t);
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        tmp[i] = v[i] + dt * k1[i];
    }
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        k2[i] = coeff * (tmp[i + 1] - 2 * tmp[i] + tmp[i - 1]) + task.f(h*i, t);
    }

}
void calc_k3(heat_task &task, const double *v, const double *k2, double *k3, double t, double *tmp)
{
    const int n = task.n;
    const double h = task.L / n;
    const double coeff = 1 / (h*h);
    const double dt = task.T / (task.m * 2.0);
    tmp[0] = task.left_condition(t);
    tmp[n] = task.right_condition(t);
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        tmp[i] = v[i] + dt * k2[i];
    }
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        k3[i] = coeff * (tmp[i + 1] - 2 * tmp[i] + tmp[i - 1]) + task.f(h*i, t);
    }
}
void calc_k4(heat_task &task, const double *v, const double *k3, double *k4, double t, double *tmp)
{
    const int n = task.n;
    const double h = task.L / n;
    const double dt = task.T / task.m;
    const double coeff = 1 / (h*h);
    tmp[0] = task.left_condition(t);
    tmp[n] = task.right_condition(t);
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        tmp[i] = v[i] + dt * k3[i];
    }
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        k4[i] = coeff * (tmp[i + 1] - 2 * tmp[i] + tmp[i - 1]) + task.f(h*i, t);
    }

}
void calc(heat_task &task, double *v, const double *k1, const double *k2, const double *k3, const double *k4, const double curTime)
{
    const int n = task.n;
    const double dt = task.T / task.m;
    const double coeff = dt / 6.0;
    #pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        v[i] += coeff * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}
void heat_equation_runge_kutta(heat_task task, double * v)//v Ц указатель на массив размера n+1, в который необходимо записать численное решение в момент t
{
    const int n = task.n;
    const double L = task.L;
    vector<double> tmp(n + 1);
    vector<double> k((n+1)*4);
    double *k1 = &k[0];
    double *k2 = &k[n+1];
    double *k3 = &k[2*n+2];
    double *k4 = &k[3*n+3];
    const int m = task.m;
    const double h = task.L / n;
    const double dt = task.T / m;
    double curTime = 0.0;
    //v(x, 0):
    for (int i = 0; i < n + 1; i++)
    {
        v[i] = task.initial_condition(i*h); //x_i = i*h
    }
    for (int tick = 0; tick < m; tick++)
    {
        curTime = dt * tick;
        v[0] = task.left_condition(curTime);
        v[n] = task.right_condition(curTime);
        calc_k1(task, v, k1, curTime);
        calc_k2(task, v, k1, k2, curTime + dt/2, tmp.data());
        calc_k3(task, v, k2, k3, curTime + dt/2, tmp.data());
        calc_k4(task, v, k3, k4, curTime + dt, tmp.data());
        calc(task, v, k1, k2, k3, k4, curTime);
    }
    v[0] = task.left_condition(curTime + dt);
    v[n] = task.right_condition(curTime + dt);
}

