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

    const int N = 3;
    double init_pt_data[N] = {1.0, 2.0, 3.0};
    pointND<N> init_pt = pointND<N>(init_pt_data);

    MPI_Datatype MPI_POINTND;
    make_pointND_MPI_Datatype(N, &MPI_POINTND);

    parareal_prob<pointND<N>> prob = parareal_prob<pointND<N>>(0.0, 1.0, init_pt, my_coarse_solve<pointND<N>>, my_fine_solve<pointND<N>>, norm_pointND<N>, MPI_POINTND);

    parareal_sol<pointND<N>> sol = parareal_sol<pointND<N>>(prob);

    parareal_sol<pointND<N>> serial_sol = parareal_sol<pointND<N>>(prob, true);

    if (my_rank == ROOT) {
        //cout << "ParaReal" << endl;
        //display_solution_csv<pointND<N>>(sol);
        //cout << "Serial" << endl;
        //display_solution_csv<pointND<N>>(serial_sol);

        cout << "ParaReal" << endl;
        display_solution_csv_2<pointND<N>>(sol);
        cout << "Serial" << endl;
        display_solution_csv_2<pointND<N>>(serial_sol);
    }

    MPI_Finalize();
    return 0;
}
