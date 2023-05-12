//#include <ctime>
#include <memory>
#include <omp.h>
#include "gtest/gtest.h"

#include "common.h"
TEST(TEST1_1, RightTest) {
  int m = 300;
  int n = 100;
  double *A = new double[m * n];
  double *B = new double[m * n];
  double *C = new double[m * n];
  double *D = new double[m * n];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i] = UniformReal(1.0, 1000.0);
      B[i] = UniformReal(1.0, 1000.0);
    }
  }

  MatAdd<double>(m, n, A, B, C);
  MatAddOmp<double>(m, n, A, B, D);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      //      EXPECT_EQ(C[i * n + j], A[i * n + j] + B[i * n + j]);
      EXPECT_DOUBLE_EQ(C[i * n + j], A[i * n + j] + B[i * n + j]);
      EXPECT_DOUBLE_EQ(D[i * n + j], A[i * n + j] + B[i * n + j]);
    }
  }
}
TEST(TEST1_1, BenchmarkTest) {
  int m = 10000;
  int n = 10000;
  int test_num =3;
  auto *a = new double[m * n];
  auto *b = new double[m * n];
  auto *c = new double[m * n];
  auto *d = new double[m * n];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      a[i] = 1.21;
      b[i] = 312.02;
    }
  }
  auto begin1 = omp_get_wtime();
  for (int i = 0; i < test_num; i++) {
    MatAdd<double>(m, n, a, b, c);
  }
  auto end1 =omp_get_wtime() - begin1;

  auto begin2 = omp_get_wtime();
  for (int i = 0; i < test_num; i++) {
    MatAddOmp<double>(m, n, a, b, d);
  }
  auto end2 = omp_get_wtime() - begin2;
  std::cout << "MatAdd time:" << end1 << " MatAddOmp time" << end2 << std::endl;
}