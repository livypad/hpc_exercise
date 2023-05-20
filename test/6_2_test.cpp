#include <memory>
#include "gtest/gtest.h"

#include "common.h"

TEST(TEST6_2, RightTest) {
  int m = 10;
  int n = 2 << m;
  auto array = new int[n];
  auto to_re_order = new int[n];

  BitReverse(n, array);
  for (int i = 0; i < n; i++) {
    EXPECT_EQ(IntReverse(i, m + 1), array[i]);
    to_re_order[i] = i;
  }
  BitReverseReOrder(n, to_re_order);
  for (int i = 0; i < n; i++) {
    EXPECT_EQ(to_re_order[i], array[i]);
  }
}