#include "common.h"
// 按位倒置排序
template <typename T>
void BitReverseReOrder(int n, T *array) {
  auto temp = new T[n];
  auto index = new int[n];
  std::memcpy(temp, array, n * sizeof(T));
  BitReverse(n, index);
  for (int i = 0; i < n; i++) {
    temp[i] = array[index[i]];
  }
  std::memcpy(array, temp, n * sizeof(T));
  delete[] temp;
  delete[] index;
}

template void BitReverseReOrder<int>(int, int*);

void BitReverse(int n, int *array) {
  if (n < 2) {
    array[0] = 0;
    return;
  }
  int n2 = 2;
  array[0] = 0;
  array[1] = 1;
  while (n2 < n) {
    for (int i = 0; i < n2; i++) {
      array[i] *= 2;
      array[i + n2] = array[i] + 1;
    }
    n2 *= 2;
  }
}

auto IntReverse(unsigned int x, int m) -> unsigned int {
  int ret = 0;
  for (int i = 0; i < m; i++) {
    if (((1 << i) & x) != 0) {
      ret |= (1 << (m - 1 - i));
    }
  }
  return ret;
}
