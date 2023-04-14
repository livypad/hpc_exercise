// 在主程序中定义了矩阵A和B,请写一个通用子程序完成矩阵相加，
// 使得对于任何m×n的子矩阵A与子矩阵B相加都能调用该子程序，
// 用C或Fortran语言（比如，matradd(m,n,,A,…,B,,C)

#include "common.h"
template <typename T>
void MatAdd(int m, int n, T *A, T *B, T *C) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      C[i * n + j] = A[i * n + j] + B[i * n + j];
    }
  }
}
template void MatAdd<double>(const int, const int, double *, double *, double *);

template <typename T>
void MatAddOmp(int m, int n, T *A, T *B, T *C) {
#pragma omp parallel for collapse(2) firstprivate(m, n, A, B, C) default(none)
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      C[i * n + j] = A[i * n + j] + B[i * n + j];
    }
  }
}
template void MatAddOmp<double>(const int, const int, double *, double *, double *);