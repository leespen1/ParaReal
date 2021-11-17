#include <mpi.h>
#include<iostream>
using std::cout; using std::endl;
#include<algorithm>
using std::copy; using std::fill_n;
#include <cmath>
using std::pow;

#define ROOT 0
/* Notes
 * Currently in the pseduocode stage, just trying to layout the structure of
 * the code to be filled out later
 *
 * Helpful resource:
 * http://fsu.digital.flvc.org/islandora/object/fsu:182428/datastream/PDF/view
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


double f(double y, double t) {
    return y;
}


double coarse_solve(double y, double t1, double t2) {
    // Single time step explicit euler
    return y + y*(t2-t1);
}


double fine_solve(double y, double t1, double t2) {
    // Mult-time step explicit euler
    int N = 25;
    double dt = (t2-t1)/N;
    for (int i=0; i < N; ++i)
        y += y*dt; 
    return y;
}


int main(int argc, char **argv){
    cout.precision(10);

    int my_rank, num_ranks;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Status status;
    /*
     * For N processes, will have N+1 points, for N fine solvers. The root
     * process will also be used as a fine solver. This works well with the
     * syntax for scatter and gather
     */
    int glo_t0 = 0;
    int glo_tf = 1;
    int num_points = num_ranks+1;

    // Define Arrays
    double *times, *points, *points_prev_rev, *points_fine_update;
    if (my_rank == ROOT) {
        cout << endl;
        times = new double[num_points];
        points = new double[num_points];
        points_prev_rev = new double[num_points];
        points_fine_update = new double[num_points];

        // Set up times
        times[0] = glo_t0;
        times[num_points-1] = glo_tf;
        double coarse_dt = double(glo_tf-glo_t0)/(num_points-1);
        for (int i=1; i < num_points-1; ++i)
            times[i] = times[i-1] + coarse_dt;

        // Set up first coarse point revision
        points[0] = 1; // Initial value
        for (int i=0; i < num_points-1; ++i)
            points[i+1] = coarse_solve(points[i], times[i], times[i+1]);
        for (int i=0; i < num_points; ++i)
            points_prev_rev[i] = points[i];

        // Check first coarse solve
        cout << "After first coarse solve:\n";
        show_points(times, points, num_points);
    }

    // Variables for this process' fine solving
    double loc_t1;
    double loc_t2;
    double loc_y1;
    double loc_y2;

    // Scatter times (these will remain constant across all iterations)
    MPI_Scatter(times, 1, MPI_DOUBLE, &loc_t1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(&times[1], 1, MPI_DOUBLE, &loc_t2, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    // Arbitrary number of revisions for now
    for (int k=0; k < 10; ++k) {
        // Scatter current revision of major points
        MPI_Scatter(points, 1, MPI_DOUBLE, &loc_y1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

        loc_y2 = fine_solve(loc_y1, loc_t1, loc_t2);

        // Gather fine points for the update
        MPI_Gather(&loc_y2, 1, MPI_DOUBLE, &points_fine_update[1], 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

        // Do the update
        if (my_rank == ROOT) {
            double next_y_prev_rev = points[1];
            double curr_y_prev_rev = points[0];
            for (int i=0; i < num_points-1; ++i) {
                next_y_prev_rev = points[i+1];
                points[i+1] = coarse_solve(points[i], times[i], times[i+1])
                              + points_fine_update[i+1]
                              - coarse_solve(points_prev_rev[i], times[i], times[i+1]);
                curr_y_prev_rev = next_y_prev_rev;
            }
            for (int i=0; i < num_points; ++i)
                points_prev_rev[i] = points[i];

            cout << "After " << k << "-th revision:\n";
            show_points(times, points, num_points);
        }
    }

    //if (my_rank == ROOT) {
    //    cout << "Final Revision\n"
    //    show_points(times, points, num_points);
    //}

    MPI_Finalize();
    return 0;
}
