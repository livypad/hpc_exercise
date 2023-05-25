// 对任意给定的x!=0和n>=0，用算法1写一个程序计算x^n。

#include "common.h"
template <typename T>
auto Xn(int n, T x) -> T {
  T result = 1;
  while (n > 0) {
    if (n % 2 == 1) {
      result *= x;
    }
    x *= x;
    n /= 2;
  }
  return result;
}

template auto Xn<double>(int , double) -> double;
template auto Xn<std::complex<double>>(int , std::complex<double>) -> std::complex<double>;