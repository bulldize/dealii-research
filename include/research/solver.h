#pragma once
#include "fields.h"
#include <deal.II/base/point.h>
#include <functional>
#include <string>
namespace research {
// Reuse the solver on the unit square with homogeneous velocity boundaries.
// Other domains/boundary types need a separate, verified mesh/BC extension.
struct Problem {
  std::string name;
  std::function<Field(const Point<2> &, double)> reference;
  std::function<Field(const Point<2> &, double)> forcing;
};
void run(const Parameters &, const Problem &, const std::string &);
} // namespace research
