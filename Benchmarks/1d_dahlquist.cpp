#include<mpi.h>
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/parareal.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/solvers.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/helpers.hpp"
#include<iostream>
using std::cout; using std::endl;
#include<cmath>
using std::pow;

#define ROOT 0
int main(int argc, char **argv) {
    cout.precision(10);

    MPI_Init(NULL,NULL);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    if (my_rank == ROOT)
        cout << endl;

    parareal_prob<double> prob = parareal_prob<double>(0.0, 1.0, 1.0, my_coarse_solve<double>, my_fine_solve<double>, std::abs, MPI_DOUBLE);

    parareal_sol<double> sol = parareal_sol<double>();
    solve_parareal<double>(prob, sol);

    parareal_sol<double> serial_sol = parareal_sol<double>();


    //parareal_sol<double> sol = parareal_sol<double>(prob);

    if (my_rank == ROOT) {
        solve_parareal_serial<double>(prob, serial_sol);
        cout << "ParaReal" << endl;
        display_solution_csv<double>(sol);
        //for (int k=0; k <= sol.num_revisions; ++k)
        //    show_points(sol.times, sol.get_pts_rev(k), sol.num_points);
        cout << "Serial" << endl;
        display_solution_csv<double>(serial_sol);
        //show_points(serial_sol.times, serial_sol.get_pts_rev(0), serial_sol.num_points);
    }
    //MPI_Datatype my_datatype = MPI::COMPLEX;
    
    MPI_Finalize();
    return 0;
}
