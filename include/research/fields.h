#pragma once
#include <array>
#include <deal.II/base/tensor.h>
#include <string>
namespace research {
using namespace dealii;
struct Parameters {
  double rho = 1, nu = 1, mu = 1, lambda = 1, dt = .001, end = .1;
  unsigned subdivisions = 16, max_newton = 15;
  double tolerance = 1e-12;
  double diffusion = -1; // negative selects dt; explicit 0 is allowed
  std::string mesh_file;
  bool contraction = false;
  double inlet_velocity = .1, length = .5;
  unsigned output_interval = 0;
  bool newtonian() const {
    return lambda == 0;
  }
  double diffusion_coefficient() const {
    return diffusion < 0 ? dt : diffusion;
  }
};
struct Field {
  std::array<double, 7> value{};
  std::array<Tensor<1, 2>, 7> gradient{};
};
inline Tensor<2, 2> deformation(const Field &s) {
  Tensor<2, 2> F;
  for (unsigned a = 0; a < 2; ++a)
    for (unsigned b = 0; b < 2; ++b)
      F[a][b] = s.value[3 + 2 * a + b];
  return F;
}
inline Tensor<2, 2> velocity_gradient(const Field &s) {
  Tensor<2, 2> g;
  for (unsigned a = 0; a < 2; ++a)
    g[a] = s.gradient[a];
  return g;
}
} // namespace research
