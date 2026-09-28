#include <deal.II/base/quadrature_lib.h>
#include <iostream>
#include <random>
#include <research/assembly.h>
#include <research/forms.h>
#include <research/manufactured.h>
#include <stdexcept>
using namespace research;
void require(bool b, const char *m) {
  if (!b)
    throw std::runtime_error(m);
}
int main() {
  try {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> d(-1, 1);
    Parameters p;
    double maxerr = 0;
    for (unsigned sample = 0; sample < 100; ++sample) {
      Field s, o, w, z, f;
      for (auto *a : {&s, &o, &w, &z, &f})
        for (unsigned k = 0; k < 7; ++k) {
          a->value[k] = d(rng);
          for (unsigned j = 0; j < 2; ++j)
            a->gradient[k][j] = d(rng);
        }
      Field plus = s, minus = s;
      const double eps = 1e-5;
      for (unsigned k = 0; k < 7; ++k) {
        plus.value[k] += eps * z.value[k];
        minus.value[k] -= eps * z.value[k];
        plus.gradient[k] += eps * z.gradient[k];
        minus.gradient[k] -= eps * z.gradient[k];
      }
      const double fd = (residual(p, plus, o, w, f) - residual(p, minus, o, w, f)) / (2 * eps),
                   an = jacobian(p, s, o, w, z);
      const double err = std::abs(fd - an) / (1 + std::abs(an));
      maxerr = std::max(maxerr, err);
      require(err < 1e-7, "Full coupled Jacobian finite-difference mismatch");
      require(std::abs(transport(o, s, s, 0, 2)) < 1e-14 &&
                  std::abs(transport(o, s, s, 3, 7)) < 1e-14,
              "Skew transport energy cancellation failed");
      const auto F = deformation(s), V = velocity_gradient(s);
      require(std::abs(scalar_product(Giesekus::stress(F), V) - scalar_product(V * F, F)) < 1e-12,
              "Elastic coupling cancellation failed");
    }
    for (unsigned sample = 0; sample < 20; ++sample) {
      Field s, o;
      for (auto *a : {&s, &o})
        for (unsigned k = 0; k < 7; ++k) {
          a->value[k] = d(rng);
          for (unsigned j = 0; j < 2; ++j)
            a->gradient[k][j] = d(rng);
        }
      PointJacobian point(s, o);
      for (unsigned c = 0; c < 7; ++c)
        for (unsigned k = 0; k < 7; ++k) {
          Shape w, z;
          w.component = c;
          z.component = k;
          w.value = d(rng);
          z.value = d(rng);
          for (unsigned j = 0; j < 2; ++j) {
            w.gradient[j] = d(rng);
            z.gradient[j] = d(rng);
          }
          double a = point.entry(p, w, z), b = jacobian(p, s, o, Field(w), Field(z));
          require(std::abs(a - b) < 1e-10 * (1 + std::abs(b)),
                  "Optimized primitive assembly differs from reference Jacobian");
        }
    }
    // Independent numerical differentiation verifies continuous manufactured PDE sources.
    for (unsigned i = 0; i < 20; ++i) {
      Point<2> x(.1 + .8 * (d(rng) + 1) / 2, .1 + .8 * (d(rng) + 1) / 2);
      double t = .05, h = 1e-5;
      auto s = exact(x, t), f = source(x, t, p);
      auto stp = exact(x, t + h), stm = exact(x, t - h);
      Tensor<1, 2> v;
      v[0] = s.value[0];
      v[1] = s.value[1];
      require(std::abs(s.gradient[0][0] + s.gradient[1][1]) < 1e-14,
              "Exact solution not divergence-free");
      for (unsigned c = 0; c < 7; ++c) {
        double time = (stp.value[c] - stm.value[c]) / (2 * h);
        double advection = 0, lap = 0, stressdiv = 0, pg = 0;
        for (unsigned j = 0; j < 2; ++j) {
          auto xp = x, xm = x;
          xp[j] += h;
          xm[j] -= h;
          auto a = exact(xp, t), b = exact(xm, t);
          double derivative = (a.value[c] - b.value[c]) / (2 * h);
          require(std::abs(derivative - s.gradient[c][j]) < 1e-7, "Manufactured gradient mismatch");
          advection += v[j] * derivative;
          lap += (a.gradient[c][j] - b.gradient[c][j]) / (2 * h);
          if (c < 2) {
            stressdiv +=
                (Giesekus::stress(deformation(a))[c][j] - Giesekus::stress(deformation(b))[c][j]) /
                (2 * h);
            if (j == c)
              pg = (a.value[2] - b.value[2]) / (2 * h);
          }
        }
        if (c < 2)
          require(std::abs(p.rho * (time + advection) - p.nu * lap + pg - p.mu * stressdiv -
                           f.value[c]) < 1e-6,
                  "Momentum source mismatch");
        if (c >= 3) {
          auto F = deformation(s);
          auto r = p.mu / (2 * p.lambda) * Giesekus::relaxation(F) - velocity_gradient(s) * F;
          require(std::abs(time + advection + r[(c - 3) / 2][(c - 3) % 2] - f.value[c]) < 1e-6,
                  "Constitutive source mismatch");
        }
      }
    }
    QWitherdenVincentSimplex<2> q(4, false);
    for (unsigned a = 0; a <= 8; ++a)
      for (unsigned b = 0; a + b <= 8; ++b) {
        double sum = 0;
        for (unsigned k = 0; k < q.size(); ++k)
          sum += q.weight(k) * std::pow(q.point(k)[0], a) * std::pow(q.point(k)[1], b);
        double exact_int = std::tgamma(a + 1) * std::tgamma(b + 1) / std::tgamma(a + b + 3);
        require(std::abs(sum - exact_int) < 1e-13, "Quadrature degree-eight exactness failed");
      }
    std::cout
        << "PASS: 100 coupled Jacobian checks (max relative error " << maxerr
        << "), skew/coupling cancellation, manufactured PDE derivatives, degree-eight quadrature\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
