// 对任意给定的x!=0和n>=0，用算法1写一个程序计算x^n。

#include "common.h"
auto Xn(int n, double x) -> double {
  double result = 1;
  while (n > 0) {
    if (n % 2 == 1) {
      result *= x;
    }
    x *= x;
    n /= 2;
  }
  return result;
}