#pragma once

#include <omp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

#define MIN(a, b) ((a) <=(b) ? (a) : (b))

class heat_task {
public:
    double T; // момент времени, в который необходимо аппроксимировать u(x, t)
    double L; // длина стержн€
    int n; // размер сетки по x
    int m; // размер сетки по t
    double initial_condition(double x); // функци€, задающа€ начальное условие
    double left_condition(double t); // функци€, задающа€ граничное условие при x = 0
    double right_condition(double t); // функци€, задающа€ граничное условие при x = L
    double f(double x, double t); // функци€, задающа€ внешнее воздействие
};
void heat_equation_runge_kutta(heat_task task, double * v)
{

}
// ‘ункци€ получает в аргументах следующие переменные :
// task Ц экземепл€р класса типа heat_task, в котором хранитс€ описание задачи.
// √арантируетс€ выполнение услови€ на сходимость метода на входных данных.
