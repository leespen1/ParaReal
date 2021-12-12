#ifndef POINTND_HPP
#define POINTND_HPP

#include<mpi.h>
#include<iostream>

#include<cmath>
using std::pow;

template <int N>
struct pointND {
    const int size = N;
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
        //double *new_data = new double[this->size];
        //pointND new_pointND = pointND<N>(new_data);
        pointND<N> new_pointND;
        for (int i = 0; i < new_pointND.size; ++i)
            new_pointND.data[i] = (*this)(i) + pt(i);
        //delete new_data;
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
        if (i >= this->size || i < 0)
            throw "Index for pointND is out of bounds";
        return this->data[i];
    }

    double& operator()(int i) {
        if (i >= this->size || i < 0)
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
        os << pt(i) << ";";
    os << pt(pt.size-1) << "]";
    return os;
}

template <int N1, int N2>
struct matrixND {
    const int size = N1*N2;
    const int num_rows = N1;
    const int num_cols = N2;
    double data[N1][N2];


    matrixND(double init_data[N1][N2]) {
        for (int i = 0; i < this->num_rows; ++i)
            for (int j = 0; j < this->num_rows; ++j)
                data[i][j] = init_data[i][j];
    }

    pointND<N1> operator*(const pointND<N2>& pt) {
        //if (this->num_cols != pt.size)
        //    throw "Matrix*Vector dimensions do not match";
        pointND<N1> new_pointND;
        int i;
        int j;
        double pt_entry;
        for (i = 0; i < this->num_rows; ++i) {
            pt_entry = 0;
            for (j=0; j < this->num_cols; ++j)
                pt_entry += (*this)(i, j) * pt(j);
            new_pointND.data[i] = pt_entry;
        }
        return new_pointND;
    }

    matrixND& operator=(const matrixND& mat) {
        this->check_size_equal(mat);
        for (int i = 0; i < this->num_rows; ++i)
            for (int j = 0; j < this->num_cols; ++j)
            this->data[i][j] = mat.data[i][j];
        return *this;
    }

    matrixND operator+(const matrixND& mat) {
        this->check_size_equal(mat);
        matrixND<N1, N2> new_matrixND;
        for (int i = 0; i < new_matrixND.size; ++i)
            new_matrixND.data[i] = (*this)(i) + mat(i);
        return new_matrixND;
    }

    matrixND operator-(const matrixND& mat) {
        this->check_size_equal(mat);
        matrixND<N1, N2> new_matrixND;
        for (int i = 0; i < new_matrixND.size; ++i)
            new_matrixND.data[i] = (*this)(i) - mat(i);
        return new_matrixND;
    }

    matrixND operator*(const double& scalar) {
        matrixND<N1, N2> new_matrixND;
        for (int i = 0; i < new_matrixND.size; ++i)
            new_matrixND.data[i] = (*this)(i)*scalar;
        return new_matrixND;
    }

    double& operator()(int i, int j) const {
        if (i >= this->num_rows || i < 0)
            throw "Row index for matrixND is out of bounds";
        if (j >= this->num_cols || j < 0)
            throw "Col index for matrixND is out of bounds";
        return this->data[i][j];
    }

    double& operator()(int i, int j) {
        if (i >= this->num_rows || i < 0)
            throw "Row index for matrixND is out of bounds";
        if (j >= this->num_cols || j < 0)
            throw "Col index for matrixND is out of bounds";
        return this->data[i][j];
    }

    double operator()(int i) const {
        if (i >= this->size || i < 0)
            throw "Index for pointND is out of bounds";
        return this->data[i];
    }

    double& operator()(int i) {
        if (i >= this->size || i < 0)
            throw "Index for pointND is out of bounds";
        return this->data[i];
    }

    void check_size_equal(const matrixND& mat) {
        if (this->size != mat.size)
            throw "Tried assignment between points of unequal size";
    }
};
#endif

