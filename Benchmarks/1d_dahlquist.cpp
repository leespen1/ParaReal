/* 1D Dahlquist Solver
 * Program for solving initial value problems of the form:
 *     u'(t, u) = lambda*u
 * , where lambda and u are real numbers.
 *
 * Can take the average time over a user-determined number of runs (given as a
 * command-line argument).
 *
 * The coarse solver used is Euler's method with 1 timestep. The fine solver is
 * also Euler's method, but using more timesteps (exact number is determined as
 * a command-line argument).
 *
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

#include<cmath>
using std::abs;

#define ROOT 0


double lambda = -1.0;

template <typename T>
T f_dahlquist(T y, double t) {
    return lambda*y;
}


int main(int argc, char *argv[]) {
    // Setup user inputs
    int num_fine_solves, num_runs, f_solves_per_rank;
    double init_val;
    double t0, tf;

    MPI_Init(&argc, &argv);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    if (argc != 7) {
        if (my_rank == ROOT)
            cout << "Usage: program_name num_fine_solve num_runs [x0 <initial value>] [t0 tf <start and end time>] [lambda <for u' = lambda*u>]\n";
        MPI_Finalize();
        return 0;
    }

    init_val = atof(argv[3]);

    t0 = atof(argv[4]);
    tf = atof(argv[5]);

    lambda = atof(argv[6]);

    num_fine_solves = atoi(argv[1]);
    num_runs = atoi(argv[2]);
    f_solves_per_rank = num_fine_solves/num_ranks; // Compute number of fine solves each process must perform

    // For computing average times across multiple runs
    const int num_points = num_ranks+1;
    double parareal_sol_durations_avg[num_points];
    fill_n(parareal_sol_durations_avg, num_points, 0.0);
    double serial_sol_duration_avg = 0;

    // Create parareal problem
    parareal_prob<double> prob = parareal_prob<double>(
        t0, tf, init_val,
        abs,
        MPI_DOUBLE
    );
    // Initialize solution variables (but do not solve yet)
    parareal_sol<double> sol;
    parareal_sol<double> serial_sol;

    // Repeatedly solve using parareal and serial
    for (int i=0; i < num_runs; ++i) {
        // Solve using parareal algorithm
        solve_parareal(
            prob, sol, 0,
            [=](double y, double t1, double t2) { // Coarse Solver, few timesteps
                return exp_euler<double>(y, t1, t2, f_dahlquist<double>, 1);
            },
            [=](double y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<double>(y, t1, t2, f_dahlquist<double>, f_solves_per_rank);
            }
        );
        // Solve using serial algorithm
        solve_parareal_serial(
            prob, serial_sol,
            [=](double y, double t1, double t2) { // Fine Solver, more timesteps
                return exp_euler<double>(y, t1, t2, f_dahlquist<double>, f_solves_per_rank);
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
        cout << "u0 = " << init_val << endl;
        cout << "t0 = " << t0 << "," << " tf = " << tf << endl; 
        cout << "u' = " << lambda << "*u\n";
        cout << "ParaReal" << endl;
        display_solution_csv_2(sol); // Apparently don't have to provide template arguments if compiler can deduce type
        cout << "Serial" << endl;
        display_solution_csv_2(serial_sol);
    }

    MPI_Finalize();
    return 0;
}

