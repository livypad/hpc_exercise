#include <memory>
#include "gtest/gtest.h"

#include "common.h"
const double PI = std::acos(-1);
TEST(TEST5_1, RightTest) {
  int n = 2 << 5;

  auto res1 = RootOfXn<double>(n);
  auto res2 = RootOfXnAc<double>(n);

  for (int i = 0; i < n; i++) {
    EXPECT_LT(fabs(res1->real() - res2->real()), 0.00001);
    EXPECT_LT(fabs(res1->imag() - res2->imag()), 0.00001);
  }
}