#pragma once
#include "fields.h"
namespace research {
struct Giesekus {
  static Tensor<2, 2> stress(const Tensor<2, 2> &F) {
    return F * transpose(F);
  }
  static Tensor<2, 2> relaxation(const Tensor<2, 2> &F) {
    return stress(F) * F - F;
  }
  static Tensor<2, 2> stress_derivative(const Tensor<2, 2> &F, const Tensor<2, 2> &H) {
    return H * transpose(F) + F * transpose(H);
  }
  static Tensor<2, 2> relaxation_derivative(const Tensor<2, 2> &F, const Tensor<2, 2> &H) {
    return stress_derivative(F, H) * F + stress(F) * H - H;
  }
};
} // namespace research
