#include<mpi.h>
#include "../Src/parareal.hpp"
#include "../Src/solvers.hpp"
#include "../Src/helpers.hpp"
#include "../Src/pointnd.hpp"
#include<iostream>
using std::cout; using std::endl;

#include<algorithm>
#include<iterator>
using std::fill_n;
using std::copy;
using std::begin; using std::end;

#define ROOT 0


double g_mat_init[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
matrixND<2,2> g_derivative_matrix = matrixND<2,2>(g_mat_init);

template <typename T>
T f_dahlquist(T y, double t) {
    return g_derivative_matrix*y;
}


int main(int argc, char *argv[]) {
    cout.precision(15); // Doubles have 15-digit precision
    int num_fine_solves, num_runs, f_solves_per_rank;
    const int N = 2;

    MPI_Init(&argc, &argv);

    if (argc != 3 && argc != 7) {
        cout << "Usage: program_name num_fine_solve num_runs [a b c d <four entries of 2x2 matrix>] " << argc << endl;
        return 0;
    }

    // Set up derivative matrix
    if (argc == 7) {
        g_mat_init[0][0] = atof(argv[3]);
        g_mat_init[0][1] = atof(argv[4]);
        g_mat_init[1][0] = atof(argv[5]);
        g_mat_init[1][1] = atof(argv[6]);
        g_derivative_matrix = matrixND<2,2>(g_mat_init);
    }

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    num_fine_solves = atoi(argv[1]);
    num_runs = atoi(argv[2]);
    f_solves_per_rank = num_fine_solves/num_ranks;

    //int num_fine_solves = 1000/num_ranks;
    double init_pt_data[N] = {1.0, 2.0};
    pointND<N> init_pt = pointND<N>(init_pt_data);

    MPI_Datatype MPI_POINTND;
    make_pointND_MPI_Datatype(N, &MPI_POINTND);

    parareal_prob<pointND<N>> prob = parareal_prob<pointND<N>>(
        0.0, 10.0, init_pt,
        norm_pointND<N>,
        MPI_POINTND
    );


    const int num_points = num_ranks+1;

    double parareal_sol_durations_avg[num_points];
    fill_n(parareal_sol_durations_avg, num_points, 0.0);
    double serial_sol_duration_avg = 0;

    parareal_sol<pointND<N>> sol;
    parareal_sol<pointND<N>> serial_sol;

    for (int i=0; i < num_runs; ++i) {

        solve_parareal(
            prob, sol, 0,
            [=](pointND<N> y, double t1, double t2) { // Coarse Solver, few timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, 1);
            },
            [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, f_solves_per_rank);
            }
        );

        solve_parareal_serial(
            prob, serial_sol,
            [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, f_solves_per_rank);
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
            //parareal_sol_durations_avg[j] /= num_runs;
            sol.sol_durations[j] = parareal_sol_durations_avg[j];
        }
        //serial_sol_duration_avg /= num_runs;
        serial_sol.sol_durations[0] = serial_sol_duration_avg;

        // Display Results
        cout << "Num Ranks: " << num_ranks << endl;
        cout << "Fine Solvers Per Rank: " << f_solves_per_rank << endl;
        cout << "u' = (" << g_derivative_matrix(0, 0) << "," << g_derivative_matrix(0, 1) << ";"
                         << g_derivative_matrix(1, 0) << "," << g_derivative_matrix(1, 1) << ")u" << endl;
        cout << "ParaReal" << endl;
        display_solution_csv_2(sol); // Apparently don't have to provide template arguments if compiler can deduce type
        cout << "Serial" << endl;
        display_solution_csv_2(serial_sol);
    }

    MPI_Finalize();
    return 0;
}

