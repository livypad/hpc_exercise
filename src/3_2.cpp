#include <cblas.h>
#include <mpi.h>
#include <cstdlib>
#include <cstring>
template <typename T>
void MatVecMul(int m, int n, T *A, T *B, T *C, T *x) {
  for (int i = 0; i < m; i++) {
    C[i] = B[i];
    for (int j = 0; j < n; j++) {
      C[i] += A[i * n + j] * x[j];
    }
  }
}
template <typename T>
void MatVecMulOmp(int m, int n, T *A, T *B, T *C, T *x) {
#pragma omp parallel for firstprivate(m, n, A, B, C, x) default(none)
  for (int i = 0; i < m; i++) {
    C[i] = B[i];
#pragma omp parallel for reduction(+ : C[i]) firstprivate(i, n, A, x) default(none)
    for (int j = 0; j < n; j++) {
      C[i] += A[i * n + j] * x[j];
    }
  }
}

template void MatVecMul<double>(int, int, double *, double *, double *, double *);
template void MatVecMulOmp<double>(int, int, double *, double *, double *, double *);

void MatVecMPI(int m, int n, double *A, double *B, double *C, double *x) {
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

  MPI_Barrier(MPI_COMM_WORLD);  // 等待进程0（root进程）的counts和displs数组设置完成
  // A，B，x需要从进程0（root进程）广播到其他进程，因为每个进程的A，B，x都是随机初始化的，一般是不同的
  MPI_Scatterv(A, counts2, displs2, MPI_DOUBLE, a_this_proc, m_this_proc * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(B, counts1, displs1, MPI_DOUBLE, b_this_proc, m_this_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MatVecMul(m_this_proc, n, a_this_proc, b_this_proc, c_this_proc, x);

  MPI_Gatherv(c_this_proc, m_this_proc, MPI_DOUBLE, C, counts1, displs1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Finalize();

  delete[] a_this_proc;
  delete[] b_this_proc;
  delete[] c_this_proc;

  if (rank != 0) {
    // 其他进程退出，不要干扰最后的判断
    exit(0);
  } else {
    delete[] displs1;
    delete[] counts1;
    delete[] displs2;
    delete[] counts2;
  }
}

void MatVecBlas(int m, int n, double *A, double *B, double *C, double *x) {
  memcpy(C, B, sizeof(double) * m);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, 1, 1.0, C, 1);
}