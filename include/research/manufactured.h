#pragma once
#include "giesekus.h"
#include <cmath>
#include <deal.II/base/point.h>
namespace research {
// Section 5.2 exact solution; continuous PDE sources exclude artificial diffusion.
inline Field exact(const Point<2> &x, double t) {
  const double a = x[0], b = x[1], e = std::exp(-t), k = 4 * std::acos(-1.0);
  const auto A = [](double z) { return z * z * (z - 1) * (z - 1); };
  const auto D = [](double z) { return 2 * z * (z - 1) * (2 * z - 1); };
  const auto DD = [](double z) { return 12 * z * z - 12 * z + 2; };
  Field s;
  s.value[0] = e * A(a) * D(b) / 2;
  s.value[1] = -e * D(a) * A(b) / 2;
  s.gradient[0][0] = e * D(a) * D(b) / 2;
  s.gradient[0][1] = e * A(a) * DD(b) / 2;
  s.gradient[1][0] = -e * DD(a) * A(b) / 2;
  s.gradient[1][1] = -e * D(a) * D(b) / 2;
  s.value[2] = e * (2 * a - 1) * (2 * b - 1);
  s.gradient[2][0] = 2 * e * (2 * b - 1);
  s.gradient[2][1] = 2 * e * (2 * a - 1);
  const double f = e * std::cos(k * a) * std::cos(k * b) / 6;
  s.value[3] = 1 + f;
  s.value[6] = 1 - f;
  s.gradient[3][0] = -k * e * std::sin(k * a) * std::cos(k * b) / 6;
  s.gradient[3][1] = -k * e * std::cos(k * a) * std::sin(k * b) / 6;
  s.gradient[6] = -s.gradient[3];
  return s;
}
inline Field source(const Point<2> &x, double t, const Parameters &p) {
  const auto s = exact(x, t);
  const double a = x[0], b = x[1], e = std::exp(-t);
  const auto A = [](double z) { return z * z * (z - 1) * (z - 1); };
  const auto D = [](double z) { return 2 * z * (z - 1) * (2 * z - 1); };
  const auto DD = [](double z) { return 12 * z * z - 12 * z + 2; };
  Tensor<1, 2> v;
  v[0] = s.value[0];
  v[1] = s.value[1];
  const std::array<double, 2> lap = {e * (DD(a) * D(b) + A(a) * (24 * b - 12)) / 2,
                                     -e * ((24 * a - 12) * A(b) + D(a) * DD(b)) / 2};
  Field f;
  for (unsigned c = 0; c < 2; ++c)
    f.value[c] = p.rho * (-s.value[c] + v * s.gradient[c]) - p.nu * lap[c] + s.gradient[2][c] -
                 p.mu * 2 * s.value[3 + 3 * c] * s.gradient[3 + 3 * c][c];
  const auto F = deformation(s);
  const auto r = p.mu / (2 * p.lambda) * Giesekus::relaxation(F) - velocity_gradient(s) * F;
  for (unsigned i = 0; i < 2; ++i)
    for (unsigned j = 0; j < 2; ++j) {
      const unsigned c = 3 + 2 * i + j;
      f.value[c] = -(s.value[c] - (i == j ? 1. : 0.)) + v * s.gradient[c] + r[i][j];
    }
  return f;
}
} // namespace research
