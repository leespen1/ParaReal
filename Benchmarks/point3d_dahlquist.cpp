#include<mpi.h>
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/parareal.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/solvers.hpp"
#include "/home/spencer/GradSchool/ParallelComputing/FinalProject/ParaReal/Src/helpers.hpp"
#include<iostream>
using std::cout; using std::endl;
#include<cmath>
using std::pow;

#define ROOT 0

struct point3D {
    double x;
    double y;
    double z;
    point3D() {} // By default do nothing
    point3D(double init_x, double init_y, double init_z) {
        x = init_x;
        y = init_y;
        z = init_z;
    }
    point3D& operator=(const point3D& pt) {
        this->x = pt.x;
        this->y = pt.y;
        this->z = pt.z;
        return *this;
    }
    point3D operator+(const point3D& pt) {
        double new_x = this->x + pt.x;
        double new_y = this->y + pt.y;
        double new_z = this->z + pt.z;
        return point3D(new_x, new_y, new_z);
    }
    point3D operator-(const point3D& pt) {
        double new_x = this->x - pt.x;
        double new_y = this->y - pt.y;
        double new_z = this->z - pt.z;
        return point3D(new_x, new_y, new_z);
    }
    point3D operator*(const double& scalar) {
        double new_x = scalar*(this->x);
        double new_y = scalar*(this->y);
        double new_z = scalar*(this->z);
        return point3D(new_x, new_y, new_z);
    }
};

double norm_point3D(point3D pt) {
    return pow(pow(pt.x, 2) + pow(pt.y, 2) + pow(pt.z, 2), 0.5);
}

std::ostream &operator<<(std::ostream &os, const point3D& pt) {
    return os << "[" << pt.x << "," << pt.y << "," << pt.z << "]";
}

int main(int argc, char **argv) {
    cout.precision(10);

    MPI_Init(NULL,NULL);

    int my_rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    if (my_rank == ROOT)
        cout << endl;

    // Define and commit the MPI_Datatype for point3D class
    MPI_Datatype MPI_POINT3D;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_POINT3D);
    MPI_Type_commit(&MPI_POINT3D);

    point3D init_pt = point3D(1.0, 2.0, 3.0);
    parareal_prob<point3D> prob = parareal_prob<point3D>(0.0, 1.0, init_pt, my_coarse_solve<point3D>, my_fine_solve<point3D>, norm_point3D, MPI_POINT3D);

    parareal_sol<point3D> sol = parareal_sol<point3D>(prob);

    parareal_sol<point3D> serial_sol = parareal_sol<point3D>(prob, true);


    if (my_rank == ROOT) {
        cout << "ParaReal" << endl;
        display_solution_csv<point3D>(sol);
        cout << "Serial" << endl;
        display_solution_csv<point3D>(serial_sol);
    }

    MPI_Finalize();
    return 0;
}
