#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST3_2, RightTest) {
  const int m = 20;
  const int n = 10;
  auto *a = new double[m * n];
  for (int i = 0; i < m * n; ++i) {
    a[i] = UniformReal(1.0, 100.0);
  }
  auto *b = new double[m];
  auto *x = new double[n];
  for (int i = 0; i < m; ++i) {
    b[i] = UniformReal(1.0, 100.0);
  }
  for (int i = 0; i < n; ++i) {
    x[i] = UniformReal(1.0, 100.0);
  }
  auto *c = new double[m];
  auto *c_omp = new double[m];
  auto *c_blas = new double[m];
  auto *c_mpi = new double[m];


  MatVecMul<double>(m, n, a, b, c, x);
  MatVecMulOmp<double>(m, n, a, b, c_omp, x);
  MatVecBlas(m, n, a, b, c_blas, x);
  MatVecMPI(m, n, a, b, c_mpi, x);

  for (int i = 0; i < m; ++i) {
    EXPECT_DOUBLE_EQ(c[i], c_omp[i]);
    EXPECT_DOUBLE_EQ(c[i], c_blas[i]);
    EXPECT_DOUBLE_EQ(c[i], c_mpi[i]);
  }
}
