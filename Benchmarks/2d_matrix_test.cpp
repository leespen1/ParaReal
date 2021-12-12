#include<mpi.h>
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/parareal.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/solvers.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/helpers.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/pointnd.hpp"
#include<iostream>
using std::cout; using std::endl;

#include<algorithm>
#include<iterator>
using std::fill_n;
using std::copy;
using std::begin; using std::end;

#define ROOT 0



template <typename T>
T f_dahlquist(T y, double t) {
    double mat_init[2][2] = {1.0, 0.0, 0.0, 1.0};
    matrixND<2,2> m = matrixND<2,2>(mat_init);
    //double lambda = -1.0;
    return m*y;
}


int main(int argc, char **argv) {
    cout.precision(15); // Doubles have 15-digit precision

    MPI_Init(NULL,NULL);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    //if (my_rank == ROOT)
    //    cout << endl;

    //int num_fine_solves = 1000/num_ranks;
    int num_fine_solves = 100/num_ranks;
    const int N = 2;
    double init_pt_data[N] = {1.0, 2.0};
    pointND<N> init_pt = pointND<N>(init_pt_data);

    MPI_Datatype MPI_POINTND;
    make_pointND_MPI_Datatype(N, &MPI_POINTND);

    parareal_prob<pointND<N>> prob = parareal_prob<pointND<N>>(
        0.0, 1.0, init_pt,
        norm_pointND<N>,
        MPI_POINTND
    );

    parareal_sol<pointND<N>> dummy_sol; // Just 

    #define NUM_RUNS 100000
    int num_points = dummy_sol.num_points;

    double parareal_sol_durations_avg[num_points];
    fill_n(parareal_sol_durations_avg, num_points, 0.0);
    double serial_sol_duration_avg = 0;

    parareal_sol<pointND<N>> sol;
    parareal_sol<pointND<N>> serial_sol;

    for (int i=0; i < NUM_RUNS; ++i) {

        solve_parareal(
            prob, sol, 0,
            [=](pointND<N> y, double t1, double t2) { // Coarse Solver, few timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, 1);
            },
            [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, num_fine_solves);
            }
        );

        solve_parareal_serial(
            prob, serial_sol,
            [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, num_fine_solves);
            }
        );

        // The following solution memory is only allocated for root
        if (my_rank == ROOT) {
            for (int j=0; j < num_points; ++j)
                parareal_sol_durations_avg[j] += sol.sol_durations[j];
            serial_sol_duration_avg += serial_sol.sol_durations[0];
        }
    }

    // Take average of times, overwrite them on the solutions (for display)
    if (my_rank == ROOT) {
        for (int j=0; j < num_points; ++j) {
            //parareal_sol_durations_avg[j] /= NUM_RUNS;
            sol.sol_durations[j] = parareal_sol_durations_avg[j];
        }
        //serial_sol_duration_avg /= NUM_RUNS;
        serial_sol.sol_durations[0] = serial_sol_duration_avg;

        // Display Results
        cout << "Num Ranks: " << num_ranks << endl;
        cout << "ParaReal" << endl;
        display_solution_csv_2(sol); // Apparently don't have to provide template arguments if compiler can deduce type
        cout << "Serial" << endl;
        display_solution_csv_2(serial_sol);
    }

    MPI_Finalize();
    return 0;
}

