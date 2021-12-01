#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/parareal.hpp"

#include<iostream>
using std::cout; using std::endl;

template <typename T>
void display_solution_csv(parareal_sol<T> sol) {
    for (int i =0; i < sol.num_points; ++i) {
        cout << sol.times[i];
        cout << ",";
    }
    cout << endl;
    for (int k=0; k <= sol.num_revisions; ++k) {
        T *pts_rev = sol.get_pts_rev(k);
        for (int i =0; i < sol.num_points; ++i) {
            cout << pts_rev[i];
            cout << ",";
        }
        cout << endl;
    }
}

/*
 * Right now this can only be ran by root, otherwise I get segfaults. I should
 * change it so that the check happens in this function, not outside.
 */
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
