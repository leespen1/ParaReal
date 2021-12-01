#ifndef PARAREAL_HPP
#define PARAREAL_HPP

/* Notes
 * Helpful resource:
 * http://fsu.digital.flvc.org/islandora/object/fsu:182428/datastream/PDF/view
 */
/* References
 * MPI Datatypes: https://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node229.htm
 */

#include <mpi.h>
#include <iostream>
using std::cout; using std::endl;

#define ROOT 0

template <typename T>
struct parareal_prob {
    double t0;
    double tf;
    T u0;
    MPI_Datatype datatype;
    T (*coarse_solve)(T, double, double);
    T (*fine_solve)(T, double, double);
    double (*norm)(T);

    parareal_prob(
        double init_t0,
        double init_tf,
        T uinit_u0,
        T (*init_coarse_solve)(T, double, double),
        T (*init_fine_solve)(T, double, double),
        double (*init_norm)(T),
        MPI_Datatype initial_datatype
    );
};


template <typename T>
struct parareal_sol {
    int num_points;
    double *times;
    int num_revisions = -1;
    T *points;

    parareal_sol();
    parareal_sol(parareal_prob<T> prob);
    T * get_pts_rev(int revision);
};


template <typename T>
void solve_parareal(
    const parareal_prob<T> &prob,
    parareal_sol<T> &sol,
    double tolerance=0.5,
    bool show_progress=true
);


template <typename T>
parareal_prob<T>::parareal_prob(
            double init_t0,
            double init_tf,
            T init_u0,
            T (*init_coarse_solve)(T, double, double),
            T (*init_fine_solve)(T, double, double),
            double (*init_norm)(T),
            MPI_Datatype initial_datatype
            )
{
        t0 = init_t0;
        tf = init_tf;
        u0 = init_u0;
        coarse_solve = init_coarse_solve;
        fine_solve = init_fine_solve;
        norm = init_norm;
        datatype = initial_datatype;
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
    MPI_Comm_size(MPI_COMM_WORLD,&num_points);
    num_points += 1;

    // Allocate memory
    times = new double[num_points];
    // Only allocate memory for points on root
    if (my_rank == ROOT)
        points = new T[num_points*num_points]; // Parareal is guarunteed to converge in at most N updates, although in this case there is no speedup
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

template <typename T>
T * parareal_sol<T>::get_pts_rev(int revision){
    return &(this->points[revision * this->num_points]);
};


/*
 * The actual parareal solver
 */
template <typename T>
void solve_parareal(
        const parareal_prob<T> &prob,
        parareal_sol<T> &sol,
        double tolerance,
        bool show_progress // Note that the default value can only be in the header file
        )
{

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Status status;

    T *points_prev_rev, *points_fine_update;

    // Set up times, helper arrays
    if (my_rank == ROOT) {
        //points_prev_rev = new T[sol.num_points];
        points_fine_update = new T[sol.num_points];

        // Set up times
        sol.times[0] = prob.t0;
        sol.times[sol.num_points-1] = prob.tf;
        double coarse_dt = double(prob.tf-prob.t0)/(sol.num_points-1);
        for (int i=1; i < sol.num_points-1; ++i)
            sol.times[i] = sol.times[i-1] + coarse_dt;

        // Set up initial value for each revision (same across all revisions)
        for (int k=0; k < sol.num_points; ++k)
            sol.get_pts_rev(k)[0] = prob.u0;

        // Set up first coarse point revision
        T *pts_init_rev = sol.get_pts_rev(0);
        for (int i=0; i < sol.num_points-1; ++i)
            pts_init_rev[i+1] = prob.coarse_solve(pts_init_rev[i], sol.times[i], sol.times[i+1]);
    }

    // Variables for this process' fine solving
    double loc_t1;
    double loc_t2;
    T loc_y1;
    T loc_y2;

    // Scatter times (these will remain constant across all iterations)
    MPI_Scatter(sol.times, 1, prob.datatype, &loc_t1, 1, prob.datatype, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(&sol.times[1], 1, prob.datatype, &loc_t2, 1, prob.datatype, ROOT, MPI_COMM_WORLD);

    bool close_enough = false;
    T *pts_prev_rev;
    T *pts_curr_rev;
    // Arbitrary number of revisions for now - should change to max num_points, but also check for final part changing by less than machine epsilon
    // Will need some sort of norm function in order to compare with machine epsilon
    for (int k = 1; k < sol.num_points; ++k) {
        sol.num_revisions = k;
        // Scatter current revision of major points
        MPI_Scatter(sol.get_pts_rev(k-1), 1, prob.datatype, &loc_y1, 1, prob.datatype, ROOT, MPI_COMM_WORLD);

        loc_y2 = prob.fine_solve(loc_y1, loc_t1, loc_t2);

        // Gather fine points for the update
        MPI_Gather(&loc_y2, 1, prob.datatype, &points_fine_update[1], 1, prob.datatype, ROOT, MPI_COMM_WORLD);

        // Do the update
        if (my_rank == ROOT) {
            pts_prev_rev = sol.get_pts_rev(k-1);
            pts_curr_rev = sol.get_pts_rev(k);

            for (int i=0; i < sol.num_points-1; ++i) {
                pts_curr_rev[i+1] = prob.coarse_solve(pts_curr_rev[i], sol.times[i], sol.times[i+1])
                              + points_fine_update[i+1]
                              - prob.coarse_solve(pts_prev_rev[i], sol.times[i], sol.times[i+1]);
            }

            double delta = prob.norm(pts_curr_rev[sol.num_points-1] - pts_prev_rev[sol.num_points-1]);
            if (delta < 0.5)
                close_enough = true;
        }
        MPI_Bcast(&close_enough, 1, MPI_CXX_BOOL, ROOT, MPI_COMM_WORLD);
        if (close_enough)
            break;
    }
}

template <typename T>
void solve_parareal_serial(
        const parareal_prob<T> &prob,
        parareal_sol<T> &sol,
        bool show_progress // Note that the default value can only be in the header file
        )
{
    sol.points[0] = prob.u0;
    sol.num_revisions = 0;
    for (int i=0; i < prob.num_points-1; ++i)
        sol.points[i+1] = prob.fine_solve(sol.points[i], sol.times[i], sol.times[i+1]);
}

#endif
