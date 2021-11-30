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
    int N = 25;
    double dt = (t2-t1)/N;
    for (int i=0; i < N; ++i)
        y += y*dt; 
    return y;
}


void show_points(double *times, double *points, int num_points) {
    cout << "t";
    for (int i=0; i < num_points; ++i)
        cout << "\t" << times[i];
    cout << endl;
    cout << "u(t)";
    for (int i=0; i < num_points; ++i)
        cout << "\t" << points[i];
    cout << endl;
}
#endif
