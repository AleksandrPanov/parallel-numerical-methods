include_directories(../RungeKutta/include)
add_executable(test_RungeKutta src/main.cpp ../RungeKutta/include/RungeKutta.h ../RungeKutta/include/print_generate.h)