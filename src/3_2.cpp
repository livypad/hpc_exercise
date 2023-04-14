#include <cblas.h>
#include <mpi.h>
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
  //  printf("rank = %d, size = %d\n", rank, size);
  int m_per_proc = m / size;
  int m_last_proc = m - m_per_proc * (size - 1);
  int m_this_proc = (rank == size - 1) ? m_last_proc : m_per_proc;
  auto *a_this_proc = new double[m_this_proc * n];
  auto *b_this_proc = new double[m_this_proc];
  auto *c_this_proc = new double[m_this_proc];

  MPI_Scatter(A, m_per_proc * n, MPI_DOUBLE, a_this_proc, m_per_proc * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(B, m_per_proc, MPI_DOUBLE, b_this_proc, m_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MatVecMul(m_this_proc, n, a_this_proc, b_this_proc, c_this_proc, x);

  MPI_Gather(c_this_proc, m_this_proc, MPI_DOUBLE, C, m_this_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Finalize();
  delete[] a_this_proc;
  delete[] b_this_proc;
  delete[] c_this_proc;
}

void MatVecBlas(int m, int n, double *A, double *B, double *C, double *x) {
  memcpy(C, B, sizeof(double) * m);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, 1, 1.0, C, 1);
}