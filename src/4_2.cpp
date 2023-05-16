#include "common.h"

template <typename T>
void TriangleMatrixSolve(T *A, T *b, T *x, int n) {
  for (int i = 0; i < n; i++) {
    T sum = 0;
    for (int j = 0; j < i; j++) {
      sum += A[i * n + j] * x[j];
    }
    x[i] = (b[i] - sum) / A[i * n + i];
  }
}
template void TriangleMatrixSolve<double>(double *, double *, double *, int);

enum { RootID = 0 };

void TriangleMatrixSolveMPI(double *A, double *b, double *x, int n) {
  MPI_Init(nullptr, nullptr);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size == 1) {
    TriangleMatrixSolve<double>(A, b, x, n);
  } else {
    int n_per_proc = (n + size - 1) / size;
    int rem = n - size * (n / size);
    int n_this_proc = rem == 0 || rank < rem ? n_per_proc : n_per_proc - 1;
    LOG_DEBUG("rank: %d, size: %d,n_this_proc: %d", rank, size, n_this_proc);

    auto displs1 = new int[size];
    auto displs2 = new int[size];
    auto counts1 = new int[size];
    auto counts2 = new int[size];

    for (int i = 0; i < size; i++) {
      displs1[i] = i == 0 ? 0 : displs1[i - 1] + counts1[i - 1];
      counts1[i] = rem == 0 || i < rem ? n : 0;

      displs2[i] = i == 0 ? 0 : displs2[i - 1] + counts2[i - 1];
      counts2[i] = rem == 0 || i < rem ? 1 : 0;
    }

    auto u = new double[n];
    auto v = new double[size - 1];

    auto a_local = new double[n * n_per_proc];
    std::memset(a_local, 0, n * n_per_proc * sizeof(double));
    auto x_local = new double[n_per_proc];

    // 按列循环发送A到各个进程
    for (int i = 0; i < n / size; i++) {
      MPI_Scatter(A + i * n * size, n, MPI_DOUBLE, a_local + i * n, n, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
    }
    if (rem != 0) {
      MPI_Scatterv(A + (n / size) * n * size, counts1, displs1, MPI_DOUBLE, a_local + (n_per_proc - 1) * n,
                   counts1[rank], MPI_DOUBLE, RootID, MPI_COMM_WORLD);
    }

    if (rank == RootID) {
      std::memcpy(u, b, n * sizeof(double));
      std::memset(v, 0, (size - 1) * sizeof(double));
    } else {
      std::memset(u, 0, n * sizeof(double));
    }

    for (int i = rank; i < size * n_per_proc; i += size) {
      if (i > 0) {
        MPI_Recv(v, size - 1, MPI_DOUBLE, ((size + i - 1) % size), MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      int local_k = i / size;
      x_local[local_k] = (u[i] + v[0]) / a_local[local_k * n + i];
      LOG_DEBUG("rank=%d,x[%d] = %lf", rank, i, x_local[local_k]);

      for (int j = 0; j < size - 2; j++) {
        v[j] = i + j + 1 < n ? v[j + 1] + u[i + j + 1] - a_local[local_k * n + i + j + 1] * x_local[local_k] : v[j + 1];
      }
      v[size - 2] = i + size - 1 < n ? u[i + size - 1] - a_local[local_k * n + i + size - 1] * x_local[local_k] : 0;

      for (int j = i + size; j < n; j++) {
        u[j] -= a_local[local_k * n + j] * x_local[local_k];
      }

      MPI_Send(v, size - 1, MPI_DOUBLE, ((i + 1) % size), 1, MPI_COMM_WORLD);
    }
    // size-1号进程的send 到RootID也需要接收，匹配上recv和send
    if (rank == 0) {
      MPI_Recv(v, size - 1, MPI_DOUBLE, (size - 1), MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //循环接收x
    for (int i = 0; i < n / size; i++) {
      MPI_Gather(x_local + i, 1, MPI_DOUBLE, x + i * size, 1, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
    }
    if (rem != 0) {
      MPI_Gatherv(x_local + n / size, rank < rem ? 1 : 0, MPI_DOUBLE, x + size * (n / size), counts2, displs2,
                  MPI_DOUBLE, RootID, MPI_COMM_WORLD);
    }
    delete[] displs1;
    delete[] displs2;
    delete[] counts1;
    delete[] counts2;
    delete[] a_local;
    delete[] u;
    delete[] v;
  }

  MPI_Finalize();

  if (rank != 0) {
    // 其他进程退出，不要干扰最后的判断
    exit(0);
  }
}