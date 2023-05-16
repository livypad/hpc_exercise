#pragma once

#include <mpi.h>
#include <cstdlib>
#include <cstring>

#include "logger.h"

auto UniformReal(double a, double b) -> double;
template <typename T>
void MatAdd(int m, int n, T *A, T *B, T *C);
template <typename T>
void MatAddOmp(int m, int n, T *A, T *B, T *C);
auto PolyResult(int n, const double *a, double x) -> double;
auto Xn(int n, double x) -> double;

template <typename T>
void MatVecMul(int m, int n, T *A, T *B, T *C, T *x);
template <typename T>
void MatVecMulOmp(int m, int n, T *A, T *B, T *C, T *X);
void MatVecBlas(int m, int n, double *A, double *B, double *C, double *x);
void MatVecMPI(int m, int n, double *A, double *B, double *C, double *x);

template <typename T>
void GaussSeidelRow(T *A, T *B, T *x, T eps, int n);
template <typename T>
void BadGaussSeidel(T *A, T *B, T *x, T eps, int n);

void RowMatrixEigen(double *A, double *B, double *x, int n);
void GaussSeidelRowMPICol(double *A, double *B, double *x, double eps, int n);

template <typename T>
void TriangleMatrixSolve(T *A, T *b, T *x, int n);
void TriangleMatrixSolveMPI(double *A, double *b, double *x, int n);