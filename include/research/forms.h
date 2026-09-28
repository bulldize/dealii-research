#pragma once
#include "giesekus.h"
namespace research {
inline double transport(const Field &old, const Field &s, const Field &w, unsigned first,
                        unsigned last) {
  Tensor<1, 2> v;
  v[0] = old.value[0];
  v[1] = old.value[1];
  double r = 0;
  for (unsigned c = first; c < last; ++c)
    r += .5 * ((v * s.gradient[c]) * w.value[c] - s.value[c] * (v * w.gradient[c]));
  return r;
}
// Paper (3.3), (3.8): old velocity is deliberately held fixed.
template <class Constitutive = Giesekus>
inline double residual(const Parameters &p, const Field &s, const Field &old, const Field &w,
                       const Field &force) {
  double r = p.rho * transport(old, s, w, 0, 2) + transport(old, s, w, 3, 7);
  for (unsigned c = 0; c < 2; ++c)
    r += p.rho * (s.value[c] - old.value[c]) * w.value[c] / p.dt +
         p.nu * (s.gradient[c] * w.gradient[c]) - force.value[c] * w.value[c];
  r -= s.value[2] * (w.gradient[0][0] + w.gradient[1][1]);
  r -= w.value[2] * (s.gradient[0][0] + s.gradient[1][1]);
  const auto F = deformation(s), G = deformation(w);
  r += (p.newtonian() ? 0. : p.mu) * scalar_product(Constitutive::stress(F), velocity_gradient(w));
  if (!p.newtonian())
    r += scalar_product(
        p.mu / (2 * p.lambda) * Constitutive::relaxation(F) - velocity_gradient(s) * F, G);
  for (unsigned c = 3; c < 7; ++c)
    r += (s.value[c] - old.value[c]) * w.value[c] / p.dt +
         p.diffusion_coefficient() * (s.gradient[c] * w.gradient[c]) - force.value[c] * w.value[c];
  return r;
}
template <class Constitutive = Giesekus>
inline double jacobian(const Parameters &p, const Field &s, const Field &old, const Field &w,
                       const Field &z) {
  double r = p.rho * transport(old, z, w, 0, 2) + transport(old, z, w, 3, 7);
  for (unsigned c = 0; c < 2; ++c)
    r += p.rho * z.value[c] * w.value[c] / p.dt + p.nu * (z.gradient[c] * w.gradient[c]);
  r -= z.value[2] * (w.gradient[0][0] + w.gradient[1][1]) +
       w.value[2] * (z.gradient[0][0] + z.gradient[1][1]);
  const auto F = deformation(s), H = deformation(z), G = deformation(w);
  r += (p.newtonian() ? 0. : p.mu) *
       scalar_product(Constitutive::stress_derivative(F, H), velocity_gradient(w));
  if (!p.newtonian())
    r += scalar_product(p.mu / (2 * p.lambda) * Constitutive::relaxation_derivative(F, H) -
                            velocity_gradient(z) * F - velocity_gradient(s) * H,
                        G);
  for (unsigned c = 3; c < 7; ++c)
    r += z.value[c] * w.value[c] / p.dt +
         p.diffusion_coefficient() * (z.gradient[c] * w.gradient[c]);
  return r;
}
} // namespace research
