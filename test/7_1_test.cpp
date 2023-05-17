#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST7_1, RightTest) {
  const int m = 20;
  const int n = 10;
  auto a = new double[m * n];
  auto a_trans = new double[m * n];
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i * n + j] = UniformReal(1.0, 10.0);
      a_trans[i + j * m] = a[i * n + j];
    }
  }
  auto *b = new double[m];
  auto *x = new double[n];
  for (int i = 0; i < m; ++i) {
    b[i] = UniformReal(1.0, 10.0);
  }
  for (int i = 0; i < n; ++i) {
    x[i] = UniformReal(1.0, 10.0);
  }
  auto *y = new double[m];
  auto *blas_y = new double[m];
  auto *mpi_y = new double[m];

  MatVecMul<double>(m, n, a, b, x, y);
  MatVecEigen(m, n, a, b, x, blas_y);
  MatVecMPICol(m, n, a_trans, b, x, mpi_y);

  for (int i = 0; i < m; ++i) {
    EXPECT_DOUBLE_EQ(y[i], blas_y[i]);
    EXPECT_DOUBLE_EQ(y[i], mpi_y[i]);
  }
}