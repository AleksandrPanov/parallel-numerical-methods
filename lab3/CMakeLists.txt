cmake_minimum_required(VERSION 2.8) #Проверка версии CMake
project(RungeKutta)			            #Название проекта
set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
add_subdirectory(RungeKutta)
add_subdirectory(test_RungeKutta)