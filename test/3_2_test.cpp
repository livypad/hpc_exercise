#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST3_2, RightTest) {
  const int m = 20;
  const int n = 10;
  auto a = new double[m * n];
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
  auto y = new double[m];
  auto omp_y = new double[m];
  auto blas_y = new double[m];
  auto mpi_y = new double[m];

  MatVecMul<double>(m, n, a, b, x, y);
  MatVecMulOmp<double>(m, n, a, b, x, omp_y);
  MatVecBlas(m, n, a, b, x, blas_y);
  MatVecMPIRow(m, n, a, b, x, mpi_y);

  for (int i = 0; i < m; ++i) {
    EXPECT_DOUBLE_EQ(y[i], omp_y[i]);
    EXPECT_DOUBLE_EQ(y[i], blas_y[i]);
    EXPECT_DOUBLE_EQ(y[i], mpi_y[i]);
  }
}
