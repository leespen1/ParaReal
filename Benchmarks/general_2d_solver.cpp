/* General 2D Parareal Solver
 * Program for solving initial value problems of the form:
 *     u'(t, u) = Au
 * , where u is a real 2D vector and A is a real 2x2 matrix (mathematically,
 * not in terms of C++ datatypes).
 *
 * Can take the average time over a user-determined number of runs (given as a
 * command-line argument).
 *
 * The coarse solver used is Euler's method with 1 timestep. The fine solver is
 * also Euler's method, but using more timesteps (exact number is determined as
 * a command-line argument).
 *
 * Default value of A is the 2x2 identity matrix, but the user choose a
 * different matrix by providing an additional 4 arguments on the command line.
 */

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
    // Setup user inputs
    int num_fine_solves, num_runs, f_solves_per_rank;
    const int N = 2;
    double init_pt_data[N];
    double t1, t2;

    MPI_Init(&argc, &argv);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    if (argc != 7 && argc != 11) {
        if (my_rank == ROOT)
            cout << "Usage: program_name num_fine_solve num_runs [x1 x2 <2 entries of initial value>] [t1 t2 <start and end time>] [a b c d <four entries of 2x2 matrix>]\n";
        MPI_Finalize();
        return 0;
    }

    init_pt_data[0] = atof(argv[3]);
    init_pt_data[1] = atof(argv[4]);
    pointND<N> init_pt = pointND<N>(init_pt_data);

    t1 = atof(argv[5]);
    t1 = atof(argv[6]);

    // Set up derivative matrix if user specified non-default
    if (argc == 11) {
        g_mat_init[0][0] = atof(argv[7]);
        g_mat_init[0][1] = atof(argv[8]);
        g_mat_init[1][0] = atof(argv[9]);
        g_mat_init[1][1] = atof(argv[10]);
    }
    g_derivative_matrix = matrixND<2,2>(g_mat_init);

    num_fine_solves = atoi(argv[1]);
    num_runs = atoi(argv[2]);
    f_solves_per_rank = num_fine_solves/num_ranks; // Compute number of fine solves each process must perform

    // Create MPI Datatype to transmit custom ND points
    MPI_Datatype MPI_POINTND;
    make_pointND_MPI_Datatype(N, &MPI_POINTND);


    // For computing average times across multiple runs
    const int num_points = num_ranks+1;
    double parareal_sol_durations_avg[num_points];
    fill_n(parareal_sol_durations_avg, num_points, 0.0);
    double serial_sol_duration_avg = 0;

    // Create parareal problem
    parareal_prob<pointND<N>> prob = parareal_prob<pointND<N>>(
        t1, t2, init_pt,
        norm_pointND<N>,
        MPI_POINTND
    );
    // Initialize solution variables (but do not solve yet)
    parareal_sol<pointND<N>> sol;
    parareal_sol<pointND<N>> serial_sol;

    // Repeatedly solve using parareal and serial
    for (int i=0; i < num_runs; ++i) {
        // Solve using parareal algorithm
        solve_parareal(
            prob, sol, 0,
            [=](pointND<N> y, double t1, double t2) { // Coarse Solver, few timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, 1);
            },
            [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, f_solves_per_rank);
            }
        );
        // Solve using serial algorithm
        solve_parareal_serial(
            prob, serial_sol,
            [=](pointND<N> y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<pointND<N>>(y, t1, t2, f_dahlquist<pointND<N>>, f_solves_per_rank);
            }
        );
        // Add times onto averages (memory only allocated on ROOT)
        if (my_rank == ROOT) {
            for (int j=0; j < num_points; ++j)
                parareal_sol_durations_avg[j] += sol.sol_durations[j];
            serial_sol_duration_avg += serial_sol.sol_durations[0];
        }
    }

    if (my_rank == ROOT) {
        // Take average of times, overwrite them on the solutions (for display)
        for (int j=0; j < num_points; ++j) {
            parareal_sol_durations_avg[j] /= num_runs;
            sol.sol_durations[j] = parareal_sol_durations_avg[j];
        }
        serial_sol_duration_avg /= num_runs;
        serial_sol.sol_durations[0] = serial_sol_duration_avg;

        // Display Results
        cout.precision(15); // Doubles have 15-digit precision
        cout << "Num Ranks: " << num_ranks << endl;
        cout << "Fine Solvers Per Rank: " << f_solves_per_rank << endl;
        cout << "u0 = (" << init_pt(0) "," << init_pt(1) << ")\n"
        cout << "t0 = " << t1 << "," << " tf = " << t2 << endl; 
        cout << "u' = (" << g_derivative_matrix(0, 0) << "," << g_derivative_matrix(0, 1) << ";"
                         << g_derivative_matrix(1, 0) << "," << g_derivative_matrix(1, 1) << ")u\n";
        cout << "ParaReal" << endl;
        display_solution_csv_2(sol); // Apparently don't have to provide template arguments if compiler can deduce type
        cout << "Serial" << endl;
        display_solution_csv_2(serial_sol);
    }

    MPI_Finalize();
    return 0;
}

