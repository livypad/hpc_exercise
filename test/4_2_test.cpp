#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST4_2, RightTest) {
  const int n = 10;
  auto *a = new double[n * n];
  auto *a_trans = new double[n * n];

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      a[i * n + j] = UniformReal(1.0, 2.0);
      a_trans[i + j * n] = a[i * n + j];
    }
  }
  auto *b = new double[n];
  auto *x = new double[n];
  auto *eigen_x = new double[n];
  auto *mpi_x = new double[n];

  for (int i = 0; i < n; ++i) {
    b[i] = UniformReal(2.0 * n, 3.0 * n);
  }

  TriangleMatrixSolve(a, b, x, n);
  RowMatrixEigen(a, b, eigen_x, n);
  TriangleMatrixSolveMPI(a_trans, b, mpi_x, n);

  for (int i = 0; i < n; ++i) {
    LOG_DEBUG("x[%d] = %lf,mpi[%d] = %lf,eigen[%d]=%lf", i, x[i], i, mpi_x[i], i, eigen_x[i]);
    EXPECT_LT(fabs(x[i] - mpi_x[i]), 0.0001);
    EXPECT_LT(fabs(x[i] - eigen_x[i]), 0.0001);
  }
}