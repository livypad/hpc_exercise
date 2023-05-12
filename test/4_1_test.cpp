#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST4_1, RightTest) {
  const int n = 3;
  auto *a = new double[n * n];
  auto *a_trans = new double[n * n];

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        a[i * n + j] = UniformReal(n * 10.0, n * 50.0);
      } else {
        a[i * n + j] = UniformReal(1.0, 10.0);
      }
      a_trans[n * j + i] = a[i * n + j];
    }
  }
  auto *b = new double[n];

  auto *mpi_x = new double[n];
  auto *bad_x = new double[n];
  auto *x = new double[n];
  auto *eigen_x = new double[n];
  for (int i = 0; i < n; ++i) {
    b[i] = UniformReal(0.0, 100.0);
  }

  RowMatrixEigen(a, b, eigen_x, n);
  for (int i = 0; i < n; ++i) {
    x[i] = eigen_x[i] + UniformReal(-50, 50);
    mpi_x[i] = x[i];
    bad_x[i] = x[i];
  }

  GaussSeidelRow<double>(a, b, x, 0.00001, n);
  for (int i = 0; i < n; ++i) {
    EXPECT_LT(fabs(x[i] - eigen_x[i]), 0.01);
  }
  
  GaussSeidelRowMPICol(a_trans, b, mpi_x, 0.000001, n);
  BadGaussSeidel<double>(a, b, bad_x, 0.000001, n);
  for (int i = 0; i < n; ++i) {
    EXPECT_DOUBLE_EQ(bad_x[i], mpi_x[i]);
    EXPECT_LT(fabs(mpi_x[i] - eigen_x[i]), 0.01);
  }
}