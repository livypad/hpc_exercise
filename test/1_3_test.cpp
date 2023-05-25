#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST1_3, RightTest) {
  const int n = 10;
  double x = UniformReal(2.0, 10.0);

  double a = Xn<double>(n, x);
  double b = pow(x, n);

  EXPECT_DOUBLE_EQ(a, b);
}
