#include <iostream>
#include <cmath>
#include <string>
#include <omp.h>
#include <ctime>
#include "RungeKutta.h"

double error(heat_task task, double * v) 
{
    double err = 0;
    for (int i = 0; i < task.n + 1; i++) 
    {
        err = max(err, task.error(v[i], i*(task.L / task.n), task.T));
    }
    return err;
}

int main() 
{
    omp_set_num_threads(4);
    heat_task task;
    if (task.m < 2 * task.T*pow((task.n / task.L), 2)) 
    {
        cout << "m should be more than " << 2 * task.T*pow((task.n / task.L), 2) << "\n";
        return 1;
    }

    double * v = new double[task.n + 1];

    int time1 = clock();
    heat_equation_runge_kutta(task, v);
    int time2 = clock() - time1;

    cout << time2 << " n = " << task.n << " , m = " << task.m << " , error = " << error(task, v) << "\n";

    return 0;
}