#pragma once
#include "forms.h"
namespace research {
// A primitive element basis occupies one component, avoiding dense zero-filled Fields.
struct Shape {
  unsigned component = 0;
  double value = 0;
  Tensor<1, 2> gradient;
  operator Field() const {
    Field f;
    f.value[component] = value;
    f.gradient[component] = gradient;
    return f;
  }
};
struct PointJacobian {
  Tensor<1, 2> old_velocity;
  Tensor<2, 2> F, grad_v;
  std::array<Tensor<2, 2>, 4> dstress, drelax;
  PointJacobian(const Field &s, const Field &old)
      : F(deformation(s)), grad_v(velocity_gradient(s)) {
    for (unsigned a = 0; a < 2; ++a)
      old_velocity[a] = old.value[a];
    for (unsigned k = 0; k < 4; ++k) {
      Tensor<2, 2> H;
      H[k / 2][k % 2] = 1;
      dstress[k] = Giesekus::stress_derivative(F, H);
      drelax[k] = Giesekus::relaxation_derivative(F, H);
    }
  }
  double entry(const Parameters &p, const Shape &w, const Shape &z) const {
    const unsigned c = w.component, d = z.component;
    double r = 0;
    if (c < 2) {
      if (d == c)
        r = p.rho / p.dt * w.value * z.value + p.nu * (w.gradient * z.gradient) +
            .5 * p.rho *
                ((old_velocity * z.gradient) * w.value - z.value * (old_velocity * w.gradient));
      if (d == 2)
        r -= z.value * w.gradient[c];
      if (d >= 3 && !p.newtonian())
        r += p.mu * z.value * (dstress[d - 3][c] * w.gradient);
    } else if (c == 2) {
      if (d < 2)
        r -= w.value * z.gradient[d];
    } else {
      const unsigned a = (c - 3) / 2, b = (c - 3) % 2;
      if (d == c)
        r = w.value * z.value / p.dt +
            .5 * ((old_velocity * z.gradient) * w.value - z.value * (old_velocity * w.gradient)) +
            p.diffusion_coefficient() * (w.gradient * z.gradient);
      if (!p.newtonian()) {
        if (d < 2 && d == a)
          for (unsigned j = 0; j < 2; ++j)
            r -= z.gradient[j] * F[j][b] * w.value;
        if (d >= 3) {
          const unsigned i = (d - 3) / 2, j = (d - 3) % 2;
          r += p.mu / (2 * p.lambda) * drelax[d - 3][a][b] * w.value * z.value;
          if (j == b)
            r -= grad_v[a][i] * w.value * z.value;
        }
      }
    }
    return r;
  }
};
} // namespace research
