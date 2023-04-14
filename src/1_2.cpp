// 设多项式Pn(x)=a[n]x^n十a[n-1]x^(n-1)十·十a[1]x十a[0],
// 请给出计算Pn(x)的方法并写出计算它的子程序。

#include "common.h"

auto PolyResult(const int n, const double *a, double x) -> double {
  double result = 0.0;
  for (int i = n; i >= 0; i--) {
    result *= x;
    result += a[i];
  }
  return result;
}