#include "common.h"
const double PI = std::acos(-1);
// n=2^m,calculate the root of x^n=1
template <typename T>
auto RootOfXn(int n) -> std::complex<T> * {
  auto result = new std::complex<T>[n];
  for (int i = 0; i < n; i++) {
    result[i] = std::polar(1.0, 2 * (i + 1) * PI / n);
  }
  return result;
}

template <typename T>
auto RootOfXnAc(int n) -> std::complex<T> * {
  std::complex<T> wn(cos(2 * PI / n), sin(2 * PI / n));
  std::complex<T> tmp(1, 0);
  auto result = new std::complex<T>[n];
  for (int i = 0; i < n; i++) {
    result[i] = tmp * wn;
    tmp *= wn;
  }
  return result;
}

template auto RootOfXn<double>(int n) -> std::complex<double> *;
template auto RootOfXnAc<double>(int n) -> std::complex<double> *;
