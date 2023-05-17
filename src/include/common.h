#pragma once

#include <mpi.h>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <memory>

#include "logger.h"

auto UniformReal(double a, double b) -> double;
template <typename T>
void MatAdd(int m, int n, T *A, T *B, T *C);
template <typename T>
void MatAddOmp(int m, int n, T *A, T *B, T *C);
auto PolyResult(int n, const double *a, double x) -> double;
auto Xn(int n, double x) -> double;

template <typename T>
void MatVecMul(int m, int n, T *A, T *B, T *x, T *y);
template <typename T>
void MatVecMulOmp(int m, int n, T *A, T *B, T *X, T *y);
void MatVecEigen(int m, int n, double *A, double *B, double *x, double *y);
void MatVecMPIRow(int m, int n, double *A, double *B, double *x, double *y);

template <typename T>
void GaussSeidelRow(T *A, T *B, T *x, T eps, int n);
template <typename T>
void BadGaussSeidel(T *A, T *B, T *x, T eps, int n);

void RowMatrixEigen(double *A, double *B, double *x, int n);
void GaussSeidelRowMPICol(double *A, double *B, double *x, double eps, int n);

template <typename T>
void TriangleMatrixSolve(T *A, T *b, T *x, int n);
void TriangleMatrixSolveMPI(double *A, double *b, double *x, int n);

template <typename T>
auto RootOfXn(int n) -> std::complex<T> *;
template <typename T>
auto RootOfXnAc(int n) -> std::complex<T> *;

void BitReverse(int n, int *array);
auto IntReverse(unsigned int x, int n) -> unsigned int;
template <typename T>
void BitReverseReOrder(int n, T *array);

void MatVecMPICol(int m, int n, double *A, double *B, double *x, double *y);