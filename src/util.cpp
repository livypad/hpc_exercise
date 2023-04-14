#include <random>
#include "common.h"
auto UniformReal(const double a, const double b) -> double {
  std::random_device rd;   // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(a, b);
  return dis(gen);
}