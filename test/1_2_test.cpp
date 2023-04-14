#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST1_2, Test1) {
  const int n = 10;

  auto *a = new double[n];
  double x = UniformReal(2.0, 10.0);
  double naive_result = 0;

  for (int i = 0; i < n; ++i) {
    a[i] = UniformReal(1.0, 100.0);
    naive_result += a[i] * pow(x, i);
  }
  auto re = PolyResult(n - 1, a, x);
  EXPECT_DOUBLE_EQ(naive_result, re);
}
