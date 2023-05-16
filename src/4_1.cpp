#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "common.h"

enum { RootID = 0 };

template <typename T>
void GaussSeidelRow(T *A, T *B, T *x, T eps, int n) {
  T *old_x = new T[n];
  T *new_x = new T[n];

  std::memcpy(new_x, x, n * sizeof(T));
  T new_norm;
  T old_norm;

  do {
    std::memcpy(old_x, new_x, n * sizeof(T));

    new_norm = 0;
    old_norm = 0;

    for (int i = 0; i < n; i++) {
      T old_sum = 0.0;
      T new_sum = 0.0;
      for (int j = 0; j < i; j++) {
        new_sum += A[n * i + j] * new_x[j];
      }
      for (int j = i + 1; j < n; j++) {
        old_sum += A[n * i + j] * old_x[j];
      }
      new_x[i] = (B[i] - old_sum - new_sum) / A[n * i + i];

      old_norm += (new_x[i] - x[i]) * (new_x[i] - x[i]);
      new_norm += (new_x[i] - old_x[i]) * (new_x[i] - old_x[i]);
    }
    old_norm = std::sqrt(old_norm);
    new_norm = std::sqrt(new_norm);
  } while (new_norm > eps * old_norm);

  std::memcpy(x, new_x, n * sizeof(T));

  delete[] old_x;
  delete[] new_x;
}

template <typename T>
void BadGaussSeidel(T *A, T *B, T *x, T eps, int n) {
  T *old_x = new T[n];
  T *new_x = new T[n];

  std::memcpy(new_x, x, n * sizeof(T));
  T new_norm;
  T old_norm;
  for (int i = 0; i < n; i++) {
    LOG_DEBUG("x[%d]=%lf", i, x[i]);
  }
  do {
    std::memcpy(old_x, new_x, n * sizeof(T));

    new_norm = 0;
    old_norm = 0;

    for (int i = 0; i < n; i++) {
      T sum = 0.0;
      for (int j = 0; j < n; j++) {
        if (j != i) {
          sum += A[n * i + j] * old_x[j];
        }
      }

      new_x[i] = (B[i] - sum) / A[n * i + i];
      LOG_DEBUG("sum: %lf,x[%d]=%lf", sum, i, x[i]);
      old_norm += (new_x[i] - x[i]) * (new_x[i] - x[i]);
      new_norm += (new_x[i] - old_x[i]) * (new_x[i] - old_x[i]);
    }
    old_norm = std::sqrt(old_norm);
    new_norm = std::sqrt(new_norm);
    LOG_DEBUG("old_norm: %lf, new_norm: %lf", old_norm, new_norm);
  } while (new_norm > eps * old_norm);

  std::memcpy(x, new_x, n * sizeof(T));

  delete[] old_x;
  delete[] new_x;
}

void GaussSeidelRowMPICol(double *A, double *B, double *x, double eps, int n) {
  MPI_Init(nullptr, nullptr);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n_per_proc = n / size;
  int n_last_proc = n - n_per_proc * (size - 1);
  int n_this_proc = (rank == size - 1) ? n_last_proc : n_per_proc;
  LOG_DEBUG("rank: %d, size: %d,n_this_proc: %d", rank, size, n_this_proc);

  auto displs1 = new int[size];
  auto displs2 = new int[size];
  auto counts1 = new int[size];
  auto counts2 = new int[size];

  for (int i = 0; i < size; i++) {
    displs1[i] = i * n_per_proc * n;
    counts1[i] = i == size - 1 ? n_last_proc * n : n_per_proc * n;

    displs2[i] = i * n_per_proc;
    counts2[i] = i == size - 1 ? n_last_proc : n_per_proc;
  }
  LOG_DEBUG("displs1: %d, counts1: %d, displs2: %d, counts2: %d", displs1[0], counts1[0], displs2[0], counts2[0]);

  auto a_local = new double[n * n_this_proc];
  auto x_local = new double[n_this_proc];

  auto x_new = new double[n_this_proc];
  auto x_old = new double[n_this_proc];

  auto x_delta = new double[n_this_proc];
  auto sum = new double[n];

  MPI_Scatterv(A, counts1, displs1, MPI_DOUBLE, a_local, n * n_this_proc, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
  MPI_Scatterv(x, counts2, displs2, MPI_DOUBLE, x_local, n_this_proc, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
  MPI_Bcast(B, n, MPI_DOUBLE, RootID, MPI_COMM_WORLD);

  double old_norm[1];
  double new_norm[1];
  std::memcpy(x_new, x_local, n_this_proc * sizeof(double));
  for (int i = 0; i < n_this_proc; i++) {
    LOG_DEBUG("x[%d]: %lf", i, x_new[i]);
  }
  do {
    std::memcpy(x_old, x_new, n_this_proc * sizeof(double));
    std::memset(sum, 0, n * sizeof(double));
    for (int j = 0; j < n_this_proc; j++) {
      for (int i = 0; i < n; i++) {
        if (i != rank * n_per_proc + j) {
          sum[i] += a_local[j * n + i] * x_old[j];
        }
      }
    }
    for (int i = 0; i < n; i++) {
      LOG_DEBUG("sum[%d]: %lf", i, sum[i]);
    }

    for (int i = 0; i < n_this_proc; i++) {
      x_new[i] = B[rank * n_per_proc + i];
    }
    //    for (int i = 0; i < n_this_proc; i++) {
    //      LOG_DEBUG("x_new[%d]: %lf", i, x_new[i]);
    //    }

    for (int i = 0; i < size; i++) {
      MPI_Scatterv(sum, counts2, displs2, MPI_DOUBLE, x_delta, n_this_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
      for (int j = 0; j < n_this_proc; j++) {
        x_new[j] -= x_delta[j];
      }
    }

    for (int i = 0; i < n_this_proc; i++) {
      x_new[i] /= a_local[i * n + i + rank * n_per_proc];
    }
    for (int i = 0; i < n_this_proc; i++) {
      LOG_DEBUG("x_new[%d]: %lf", i, x_new[i]);
    }

    old_norm[0] = 0;
    new_norm[0] = 0;
    for (int i = 0; i < n_this_proc; i++) {
      old_norm[0] += (x_new[i] - x_local[i]) * (x_new[i] - x_local[i]);
      new_norm[0] += (x_new[i] - x_old[i]) * (x_new[i] - x_old[i]);
    }
    MPI_Allreduce(old_norm, old_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(new_norm, new_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    old_norm[0] = std::sqrt(old_norm[0]);
    new_norm[0] = std::sqrt(new_norm[0]);
    LOG_DEBUG("old_norm: %lf, new_norm: %lf", old_norm[0], new_norm[0]);
  } while (new_norm[0] > eps * old_norm[0]);
  for (int i = 0; i < n_this_proc; i++) {
    LOG_DEBUG("x_new[%d]: %lf", i, x_new[i]);
  }

  MPI_Gatherv(x_new, n_this_proc, MPI_DOUBLE, x, counts2, displs2, MPI_DOUBLE, RootID, MPI_COMM_WORLD);

  MPI_Finalize();

  if (rank != 0) {
    // 其他进程退出，不要干扰最后的判断
    exit(0);
  } else {
    delete[] a_local;
    delete[] x_local;
    delete[] x_new;
    delete[] x_old;
    delete[] displs1;
    delete[] counts1;
    delete[] displs2;
    delete[] counts2;
  }
}

template void GaussSeidelRow<double>(double *, double *, double *, double, int);
template void BadGaussSeidel<double>(double *, double *, double *, double, int);

void RowMatrixEigen(double *A, double *B, double *x, int n) {
  auto mat_a = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(A, n, n);
  auto vec_b = Eigen::Map<Eigen::Vector<double, Eigen::Dynamic>>(B, n);

  Eigen::Vector<double, Eigen::Dynamic> result = mat_a.inverse() * vec_b;

  std::memcpy(x, result.data(), n * sizeof(double));
}