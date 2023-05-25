#include <Eigen/Core>
#include <Eigen/Dense>
#include "include/common.h"

template <typename T>
void MatVecMul(int m, int n, T *A, T *B, T *x, T *y) {
  for (int i = 0; i < m; i++) {
    y[i] = B[i];
    for (int j = 0; j < n; j++) {
      y[i] += A[i * n + j] * x[j];
    }
  }
}
template <typename T>
void MatVecMulOmp(int m, int n, T *A, T *B, T *x, T *y) {
#pragma omp parallel for firstprivate(m, n, A, B, y, x) default(none)
  for (int i = 0; i < m; i++) {
    y[i] = B[i];
#pragma omp parallel for reduction(+ : y[i]) firstprivate(i, n, A, x) default(none)
    for (int j = 0; j < n; j++) {
      y[i] += A[i * n + j] * x[j];
    }
  }
}

template void MatVecMul<double>(int, int, double *, double *, double *, double *);
template void MatVecMulOmp<double>(int, int, double *, double *, double *, double *);

enum { RootId = 0 };

void MatVecMPIRow(int m, int n, double *A, double *B, double *x, double *y) {
  MPI_Init(nullptr, nullptr);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int m_per_proc = m / size;
  int m_last_proc = m - m_per_proc * (size - 1);
  int m_this_proc = (rank == size - 1) ? m_last_proc : m_per_proc;

  int *displs1;
  int *displs2;
  int *counts1;
  int *counts2;

  if (rank == 0) {
    displs1 = new int[size];
    counts1 = new int[size];
    displs2 = new int[size];
    counts2 = new int[size];
    for (int i = 0; i < size; i++) {
      displs1[i] = i * m_per_proc;
      counts1[i] = i == size - 1 ? m_last_proc : m_per_proc;

      displs2[i] = i * m_per_proc * n;
      counts2[i] = i == size - 1 ? m_last_proc * n : m_per_proc * n;
    }
  }

  auto *a_this_proc = new double[m_this_proc * n];
  auto *b_this_proc = new double[m_this_proc];
  auto *c_this_proc = new double[m_this_proc];

  // 等待进程0（root进程）的counts和displs数组设置完成
  // A，B，x需要从进程0（root进程）广播到其他进程，因为每个进程的A，B，x都是随机初始化的，一般是不同的
  MPI_Scatterv(A, counts2, displs2, MPI_DOUBLE, a_this_proc, m_this_proc * n, MPI_DOUBLE, RootId, MPI_COMM_WORLD);
  MPI_Scatterv(B, counts1, displs1, MPI_DOUBLE, b_this_proc, m_this_proc, MPI_DOUBLE, RootId, MPI_COMM_WORLD);
  MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MatVecMul(m_this_proc, n, a_this_proc, b_this_proc, x, c_this_proc);

  MPI_Gatherv(c_this_proc, m_this_proc, MPI_DOUBLE, y, counts1, displs1, MPI_DOUBLE, RootId, MPI_COMM_WORLD);
  MPI_Finalize();

  delete[] a_this_proc;
  delete[] b_this_proc;
  delete[] c_this_proc;

  if (rank != RootId) {
    // 其他进程退出，不要干扰最后的判断
    exit(0);
  } else {
    delete[] displs1;
    delete[] counts1;
    delete[] displs2;
    delete[] counts2;
  }
}

void MatVecEigen(int m, int n, double *A, double *B, double *x, double *y) {
  auto mat_a = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(A, m, n);
  auto vec_x = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(x, n);
  auto vec_b = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(B, m);
  Eigen::Matrix<double, Eigen::Dynamic, 1> result = mat_a * vec_x + vec_b;
  std::memcpy(y, result.data(), m * sizeof(double));
}