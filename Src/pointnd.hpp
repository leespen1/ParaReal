#ifndef POINTND_HPP
#define POINTND_HPP

#include<mpi.h>
#include<iostream>

#include<cmath>
using std::pow;

template <int N>
struct pointND {
    int size = N;
    double data[N];

    pointND() {};
    pointND(double init_data[N]) {
        for (int i = 0; i < this->size; ++i)
            data[i] = init_data[i];
    }

    pointND& operator=(const pointND& pt) {
        this->check_size_equal(pt);
        for (int i = 0; i < this->size; ++i)
            this->data[i] = pt.data[i];
        return *this;
    }

    pointND operator+(const pointND& pt) {
        this->check_size_equal(pt);
        double *new_data = new double[this->size];
        pointND new_pointND = pointND<N>(new_data);
        for (int i = 0; i < new_pointND.size; ++i)
            new_pointND.data[i] = (*this)(i) + pt(i);
        return new_pointND;
    }

    pointND operator-(const pointND& pt) {
        this->check_size_equal(pt);
        double *new_data = new double[this->size];
        pointND new_pointND = pointND<N>(new_data);
        for (int i = 0; i < new_pointND.size; ++i)
            new_pointND.data[i] = (*this)(i) - pt(i);
        return new_pointND;
    }

    pointND operator*(const double& scalar) {
        double *new_data = new double[this->size];
        pointND new_pointND = pointND<N>(new_data);
        for (int i = 0; i < new_pointND.size; ++i)
            new_pointND.data[i] = (*this)(i)*scalar;
        return new_pointND;
    }

    // Not sure why I need both overloads of the '()' operator
    double operator()(int i) const {
        if (i > this->size || i < 0)
            throw "Index for pointND is out of bounds";
        return this->data[i];
    }

    double& operator()(int i) {
        if (i > this->size || i < 0)
            throw "Index for pointND is out of bounds";
        return this->data[i];
    }

    void check_size_equal(const pointND& pt) {
        if (this->size != pt.size)
            throw "Tried assignment between points of unequal size";
    }

};


void make_pointND_MPI_Datatype(int N, MPI_Datatype *data_type_ptr) {
    // Define and commit the MPI_Datatype for pointND class
    int block_lengths[2] = {1, N};
    MPI_Aint displacements[2] = {0, sizeof(double)}; // struct is padded, so need to leave the integer the space of a double
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Type_create_struct(2, block_lengths, displacements, types, data_type_ptr);
    MPI_Type_commit(data_type_ptr);
}


// The L-2 norm
template <int N>
double norm_pointND(pointND<N> pt) {
    double norm = 0;
    for (int i=0; i < pt.size; ++i)
        norm += pow(pt(i), 2);
    norm = pow(norm, 0.5);
    return norm;
}


template <int N>
std::ostream &operator<<(std::ostream &os, const pointND<N>& pt) {
    os << "[";
    for (int i=0; i < pt.size-1; ++i)
        os << pt(i) << ",";
    os << pt(pt.size-1) << "]";
    return os;
}

#endif
