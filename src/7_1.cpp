#include "common.h"

enum { RootId = 0 };

void MatVecMPICol(int m, int n, double *A, double *B, double *x, double *y) {
  MPI_Init(nullptr, nullptr);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n_per_proc = n / size;
  int n_last_proc = n - n_per_proc * (size - 1);
  int n_this_proc = (rank == size - 1) ? n_last_proc : n_per_proc;

  int *displs1;
  int *displs2;
  int *counts1;
  int *counts2;

  if (rank == RootId) {
    displs1 = new int[size];
    counts1 = new int[size];
    displs2 = new int[size];
    counts2 = new int[size];
    for (int i = 0; i < size; i++) {
      displs1[i] = i * n_per_proc;
      counts1[i] = i == size - 1 ? n_last_proc : n_per_proc;

      displs2[i] = i * n_per_proc * m;
      counts2[i] = i == size - 1 ? n_last_proc * m : n_per_proc * m;
    }
    std::memset(y, 0, m * sizeof(double));
  }

  auto *a_local = new double[n_this_proc * m];
  auto *x_local = new double[n_this_proc];
  auto *y_local = new double[m];
  std::memset(y_local, 0, m * sizeof(double));

  // A，x需要从进程0（root进程）广播到其他进程，因为每个进程的A，x都是随机初始化的，一般是不同的
  MPI_Scatterv(A, counts2, displs2, MPI_DOUBLE, a_local, n_this_proc * m, MPI_DOUBLE, RootId, MPI_COMM_WORLD);
  MPI_Scatterv(x, counts1, displs1, MPI_DOUBLE, x_local, n_this_proc, MPI_DOUBLE, RootId, MPI_COMM_WORLD);

  for (int i = 0; i < n_this_proc; i++) {
    for (int j = 0; j < m; j++) {
      y_local[j] += a_local[i * m + j] * x_local[i];
    }
  }

  MPI_Reduce(y_local, y, m, MPI_DOUBLE, MPI_SUM, RootId, MPI_COMM_WORLD);

  MPI_Finalize();

  delete[] a_local;
  delete[] x_local;
  delete[] y_local;
  delete[] displs1;
  delete[] counts1;
  delete[] displs2;
  delete[] counts2;

  if (rank != 0) {
    // 其他进程退出，不要干扰最后的判断
    exit(0);
  } else {
    // 最后在Root进程，只加一遍的B[]
    for (int i = 0; i < m; i++) {
      y[i] += B[i];
    }
  }
}
