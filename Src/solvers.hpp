#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include<iostream>
using std::cout; using std::endl;

template <typename T>
T my_coarse_solve(T y, double t1, double t2) {
    // Single time step explicit euler
    return y + y*(t2-t1);
}


template <typename T>
T my_fine_solve(T y, double t1, double t2) {
    // Mult-time step explicit euler
    int N = 100;
    double dt = (t2-t1)/N;
    for (int i=0; i < N; ++i)
        y = y + y*dt; // Not using += means I have one less operator to overload
    return y;
}

template <typename T>
T exp_euler(T &y, double &t1, double &t2, int N=1) {
    // Mult-time step explicit euler
    double dt = (t2-t1)/N;
    for (int i=0; i < N; ++i)
        y = y + y*dt; // Not using += means I have one less operator to overload
    return y;
}

#endif
