#include <mpi.h>
#include<iostream>
using std::cout; using std::endl;
#include<algorithm>
using std::copy; using std::fill_n;
#include <cmath>
using std::pow;
#include<utility>
using std::pair;
#include<tuple>
using std::tuple;

#define ROOT 0
/* Notes
 * Helpful resource:
 * http://fsu.digital.flvc.org/islandora/object/fsu:182428/datastream/PDF/view
 */

template <typename T>
struct parareal_prob {
    double t0;
    double tf;
    T (*coarse_solve)(T, double, double);
    T (*fine_solve)(T, double, double);
    parareal_prob(
            double init_t0,
            double init_tf,
            T (*init_coarse_solve)(T, double, double),
            T (*init_fine_solve)(T, double, double)
            ){
        t0 = init_t0;
        tf = init_tf;
        coarse_solve = init_coarse_solve;
        fine_solve = init_fine_solve;
    }
};


template <typename T>
struct parareal_sol {
    int num_points;
    double *times;
    T *points;
    parareal_sol();
    parareal_sol(parareal_prob<T> prob);
    // Do I need this?
    /*
    ~parareal_sol() {
        delete[] times;
        delete[] points; 
    }
    */
    /*
    parareal_sol<T>& operator=(parareal_sol<T> other) {
        num_points = other.num_points;
        times = other.times;
        points = other.points;
        return *this;
    }
    */
};


/*
 * Default constructor. I don't think I need anything else.
 *
 * For N processes, will have N+1 points, for N fine solvers. The root
 * process will also be used as a fine solver. This works well with the
 * syntax for scatter and gather. Fine solvers will have num_points=0, since
 * they do not need to allocate memory for the major points.
 */
template <typename T>
parareal_sol<T>::parareal_sol(){
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if (my_rank != ROOT)
        num_points = 0;
    else {
        MPI_Comm_size(MPI_COMM_WORLD,&num_points);
        num_points += 1;
    }
    times = new double[num_points];
    points = new T[num_points];
};

/*
 * A way to solve a parareal problem right off the bat
 */
template <typename T>
parareal_sol<T>::parareal_sol(parareal_prob<T> prob){
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if (my_rank != ROOT)
        num_points = 0;
    else {
        MPI_Comm_size(MPI_COMM_WORLD,&num_points);
        num_points += 1;
    }
    times = new double[num_points];
    points = new T[num_points];
    solve_parareal(prob, *this);
};


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


double my_coarse_solve(double y, double t1, double t2) {
    // Single time step explicit euler
    return y + y*(t2-t1);
}


double my_fine_solve(double y, double t1, double t2) {
    // Mult-time step explicit euler
    int N = 25;
    double dt = (t2-t1)/N;
    for (int i=0; i < N; ++i)
        y += y*dt; 
    return y;
}

/*
 * The actual parareal solver
 */
template <typename T>
void solve_parareal(
        const parareal_prob<T> &prob,
        parareal_sol<T> &sol,
        bool show_progress=false
        ) {

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Status status;


    T *points_prev_rev, *points_fine_update;

    // Set up times, helper arrays
    if (my_rank == ROOT) {
        points_prev_rev = new T[sol.num_points];
        points_fine_update = new T[sol.num_points];

        // Set up times
        sol.times[0] = prob.t0;
        sol.times[sol.num_points-1] = prob.tf;
        double coarse_dt = double(prob.tf-prob.t0)/(sol.num_points-1);
        for (int i=1; i < sol.num_points-1; ++i)
            sol.times[i] = sol.times[i-1] + coarse_dt;

        // Set up first coarse point revision
        sol.points[0] = 1; // Initial value
        for (int i=0; i < sol.num_points-1; ++i)
            sol.points[i+1] = prob.coarse_solve(sol.points[i], sol.times[i], sol.times[i+1]);
        for (int i=0; i < sol.num_points; ++i)
            points_prev_rev[i] = sol.points[i];

        // Check first coarse solve
        if (show_progress) {
            cout << "After first coarse solve:\n";
            show_points(sol.times, sol.points, sol.num_points);
        }
    }

    // Variables for this process' fine solving
    double loc_t1;
    double loc_t2;
    double loc_y1;
    double loc_y2;

    // Scatter times (these will remain constant across all iterations)
    MPI_Scatter(sol.times, 1, MPI_DOUBLE, &loc_t1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(&sol.times[1], 1, MPI_DOUBLE, &loc_t2, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    // Arbitrary number of revisions for now
    for (int k=0; k < 10; ++k) {
        // Scatter current revision of major points
        MPI_Scatter(sol.points, 1, MPI_DOUBLE, &loc_y1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

        loc_y2 = prob.fine_solve(loc_y1, loc_t1, loc_t2);

        // Gather fine points for the update
        MPI_Gather(&loc_y2, 1, MPI_DOUBLE, &points_fine_update[1], 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

        // Do the update
        if (my_rank == ROOT) {
            double next_y_prev_rev = sol.points[1];
            double curr_y_prev_rev = sol.points[0];
            for (int i=0; i < sol.num_points-1; ++i) {
                next_y_prev_rev = sol.points[i+1];
                sol.points[i+1] = prob.coarse_solve(sol.points[i], sol.times[i], sol.times[i+1])
                              + points_fine_update[i+1]
                              - prob.coarse_solve(points_prev_rev[i], sol.times[i], sol.times[i+1]);
                curr_y_prev_rev = next_y_prev_rev;
            }
            for (int i=0; i < sol.num_points; ++i)
                points_prev_rev[i] = sol.points[i];

            if (show_progress) {
                cout << "After " << k << "-th revision:\n";
                show_points(sol.times, sol.points, sol.num_points);
            }
        }
    }
}


int main(int argc, char **argv) {
    cout.precision(10);

    MPI_Init(NULL,NULL);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    if (my_rank == ROOT)
        cout << endl;

    parareal_prob<double> prob = parareal_prob<double>(0.0, 1.0, my_coarse_solve, my_fine_solve);
    //parareal_sol<double> sol = parareal_sol<double>();
    //solve_parareal<double>(prob, sol);
    parareal_sol<double> sol = parareal_sol<double>(prob);

    if (my_rank == ROOT) {
        cout << "Number of Points: " << sol.num_points << endl;
        cout << "Final Revision:\n";
        show_points(sol.times, sol.points, sol.num_points);
    }

    MPI_Finalize();
    return 0;
}
