#include<mpi.h>
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/parareal.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/solvers.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/helpers.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/pointnd.hpp"
#include<iostream>
using std::cout; using std::endl;

#include<algorithm>
#include<iterator>
using std::copy;
using std::begin; using std::end;

#define ROOT 0


int main(int argc, char **argv) {
    cout.precision(15); // Doubles have 15-digit precision

    MPI_Init(NULL,NULL);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    //if (my_rank == ROOT)
    //    cout << endl;

    int num_fine_solves = 1000/num_ranks;
    const int N = 3;
    double init_pt_data[N] = {1.0, 2.0, 3.0};
    pointND<N> init_pt = pointND<N>(init_pt_data);

    MPI_Datatype MPI_POINTND;
    make_pointND_MPI_Datatype(N, &MPI_POINTND);

    parareal_prob<pointND<N>> prob = parareal_prob<pointND<N>>(
        0.0, 1.0, init_pt,
        norm_pointND<N>,
        MPI_POINTND
    );

    parareal_sol<pointND<N>> sol = parareal_sol<pointND<N>>();
    parareal_sol<pointND<N>> serial_sol = parareal_sol<pointND<N>>();

    solve_parareal(
        prob, sol, 0,
        [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
            return exp_euler<pointND<N>>(y, t1, t2, 1);
        },
        [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
            return exp_euler<pointND<N>>(y, t1, t2, num_fine_solves);
        }
    );

    solve_parareal_serial(
        prob, serial_sol,
        [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
            return exp_euler<pointND<N>>(y, t1, t2, num_fine_solves);
        }
    );

    if (my_rank == ROOT) {
        //cout << "ParaReal" << endl;
        //display_solution_csv<pointND<N>>(sol);
        //cout << "Serial" << endl;
        //display_solution_csv<pointND<N>>(serial_sol);

        cout << "ParaReal" << endl;
        display_solution_csv_2(sol); // Apparently don't have to provide template arguments if compiler can deduce type
        cout << "Serial" << endl;
        display_solution_csv_2(serial_sol);
    }

    MPI_Finalize();
    return 0;
}


