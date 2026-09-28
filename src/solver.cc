#include <chrono>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <research/assembly.h>
#include <research/forms.h>
#include <research/solver.h>
#include <sstream>
namespace research {
namespace {
// Friedrichs-Keller grid: two consistently oriented triangles per square.
void square_mesh(Triangulation<2> &tria, unsigned n) {
  std::vector<Point<2>> vertices;
  for (unsigned j = 0; j <= n; ++j)
    for (unsigned i = 0; i <= n; ++i)
      vertices.emplace_back(double(i) / n, double(j) / n);
  std::vector<CellData<2>> cells;
  for (unsigned j = 0; j < n; ++j)
    for (unsigned i = 0; i < n; ++i) {
      const unsigned a = j * (n + 1) + i, b = a + 1, c = a + n + 1, d = c + 1;
      CellData<2> t(3);
      t.vertices = {a, b, d};
      cells.push_back(t);
      t.vertices = {a, d, c};
      cells.push_back(t);
    }
  tria.create_triangulation(vertices, cells, SubCellData());
}
class Inlet : public Function<2> {
public:
  Inlet(double v, double l) : Function<2>(7), v(v), l(l) {}
  double value(const Point<2> &x, unsigned c) const override {
    if (c == 0)
      return v * (1 - std::pow(x[1] / (4 * l), 2));
    return (c == 3 || c == 6) ? 1. : 0.;
  }

private:
  double v, l;
};
struct FaceCache {
  std::vector<std::vector<Shape>> shapes;
  std::vector<double> weights;
  std::vector<Tensor<1, 2>> normals;
  std::vector<types::global_dof_index> indices;
};
struct BasisCache {
  std::vector<std::vector<Shape>> shapes;
  std::vector<double> weights;
  std::vector<Point<2>> points;
  std::vector<types::global_dof_index> indices;
};
Field evaluate(const BasisCache &c, unsigned q, const Vector<double> &x) {
  Field r;
  for (unsigned i = 0; i < c.indices.size(); ++i) {
    const double a = x[c.indices[i]];
    const auto &basis = c.shapes[q][i];
    r.value[basis.component] += a * basis.value;
    r.gradient[basis.component] += a * basis.gradient;
  }
  return r;
}
} // namespace
void run(const Parameters &p, const Problem &problem, const std::string &out) {
  const auto start = std::chrono::steady_clock::now();
  const unsigned steps = std::lround(p.end / p.dt);
  AssertThrow(steps > 0 && std::abs(steps * p.dt - p.end) < 1e-12,
              ExcMessage("end_time must be an integer multiple of dt"));
  AssertThrow(p.newtonian() || p.dt < p.lambda / p.mu,
              ExcMessage("Use dt < lambda/mu for the paper's stability condition"));
  std::filesystem::create_directories(out);
  AssertThrow(!std::filesystem::exists(out + "/history.csv"),
              ExcMessage("Output exists; choose a fresh directory to preserve results"));
  Triangulation<2> tria;
  if (p.contraction) {
    AssertThrow(!p.mesh_file.empty(), ExcMessage("contraction requires mesh_file"));
    GridIn<2> reader;
    reader.attach_triangulation(tria);
    std::ifstream input(p.mesh_file);
    AssertThrow(input.good(), ExcMessage("Cannot open mesh"));
    reader.read_msh(input);
    for (const auto &cell : tria.active_cell_iterators())
      for (const auto f : cell->face_indices())
        if (cell->face(f)->at_boundary()) {
          const auto x = cell->face(f)->center();
          cell->face(f)->set_boundary_id(std::abs(x[0] + 20 * p.length) < 1e-8
                                             ? 1
                                             : (std::abs(x[0] - 40 * p.length) < 1e-8 ? 2 : 0));
        }
  } else
    square_mesh(tria, p.subdivisions);
  FESystem<2> fe(FE_SimplexP<2>(2), 2, FE_SimplexP<2>(1), 5);
  MappingFE<2> mapping{FE_SimplexP<2>(1)};
  // Degree-eight exact symmetric simplex quadrature, same degree as Section 5.2.
  QWitherdenVincentSimplex<2> quadrature(4, false);
  DoFHandler<2> dh(tria);
  dh.distribute_dofs(fe);
  AffineConstraints<double> constraints, solution_constraints;
  ComponentMask velocity_mask(7, false);
  velocity_mask.set(0, true);
  velocity_mask.set(1, true);
  std::map<types::global_dof_index, double> boundary;
  VectorTools::interpolate_boundary_values(mapping, dh, 0, Functions::ZeroFunction<2>(7), boundary,
                                           velocity_mask);
  if (p.contraction) {
    ComponentMask inlet_mask(7, true);
    inlet_mask.set(2, false);
    VectorTools::interpolate_boundary_values(mapping, dh, 1, Inlet(p.inlet_velocity, p.length),
                                             boundary, inlet_mask);
  }
  const auto pressure_dofs =
      DoFTools::extract_dofs(dh, fe.component_mask(FEValuesExtractors::Scalar(2)));
  if (!p.contraction)
    boundary[*pressure_dofs.begin()] = 0;
  if (p.newtonian())
    for (unsigned c = 3; c < 7; ++c) {
      const auto indices =
          DoFTools::extract_dofs(dh, fe.component_mask(FEValuesExtractors::Scalar(c)));
      for (const auto i : indices)
        boundary[i] = (c == 3 || c == 6) ? 1. : 0.;
    }
  for (const auto &[i, value] : boundary) {
    constraints.add_line(i);
    solution_constraints.add_line(i);
    solution_constraints.set_inhomogeneity(i, value);
  }
  constraints.close();
  solution_constraints.close();
  DynamicSparsityPattern dsp(dh.n_dofs());
  DoFTools::make_sparsity_pattern(dh, dsp, constraints, false);
  SparsityPattern sparsity;
  sparsity.copy_from(dsp);
  SparseMatrix<double> matrix(sparsity);
  Vector<double> solution(dh.n_dofs()), old(dh.n_dofs()), rhs(dh.n_dofs()), increment(dh.n_dofs());
  FEValues<2> fv(mapping, fe, quadrature,
                 update_values | update_gradients | update_JxW_values | update_quadrature_points);
  const unsigned nd = fe.n_dofs_per_cell();
  std::vector<BasisCache> cache;
  for (const auto &cell : dh.active_cell_iterators()) {
    fv.reinit(cell);
    BasisCache c;
    c.indices.resize(nd);
    cell->get_dof_indices(c.indices);
    c.weights = fv.get_JxW_values();
    c.points = fv.get_quadrature_points();
    c.shapes.resize(quadrature.size(), std::vector<Shape>(nd));
    for (unsigned q = 0; q < quadrature.size(); ++q)
      for (unsigned i = 0; i < nd; ++i) {
        const unsigned k = fe.system_to_component_index(i).first;
        c.shapes[q][i].component = k;
        c.shapes[q][i].value = fv.shape_value(i, q);
        c.shapes[q][i].gradient = fv.shape_grad(i, q);
      }
    cache.push_back(std::move(c));
  }
  std::vector<FaceCache> faces;
  if (p.contraction) {
    QGauss<1> fq(5);
    FEFaceValues<2> ff(mapping, fe, fq, update_values | update_JxW_values | update_normal_vectors);
    for (const auto &cell : dh.active_cell_iterators())
      for (const auto f : cell->face_indices())
        if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == 2) {
          ff.reinit(cell, f);
          FaceCache c;
          c.indices.resize(nd);
          cell->get_dof_indices(c.indices);
          c.weights = ff.get_JxW_values();
          c.normals = ff.get_normal_vectors();
          c.shapes.resize(fq.size(), std::vector<Shape>(nd));
          for (unsigned q = 0; q < fq.size(); ++q)
            for (unsigned i = 0; i < nd; ++i)
              {c.shapes[q][i].component=fe.system_to_component_index(i).first; c.shapes[q][i].value=ff.shape_value(i,q);}
          faces.push_back(std::move(c));
        }
  }
  auto solve_linear = [&](Vector<double> &x) {
    SparseDirectUMFPACK lu;
    lu.initialize(matrix);
    lu.vmult(x, rhs);
    constraints.distribute(x);
  };
  // L2 projection into the discretely divergence-free velocity subspace.
  // Pressure is a projection multiplier; it is not time differentiated.
  if (!p.contraction) {
    matrix = 0;
    rhs = 0;
    for (const auto &c : cache) {
      FullMatrix<double> local(nd, nd);
      Vector<double> b(nd);
      for (unsigned q = 0; q < quadrature.size(); ++q) {
        const auto ex = problem.reference(c.points[q], 0);
        for (unsigned i = 0; i < nd; ++i) {
          const Field w = c.shapes[q][i];
          for (unsigned k = 0; k < 7; ++k)
            if (k != 2)
              b[i] += ex.value[k] * w.value[k] * c.weights[q];
          for (unsigned j = 0; j < nd; ++j) {
            const Field z = c.shapes[q][j];
            double a = 0;
            for (unsigned k = 0; k < 7; ++k)
              if (k != 2)
                a += z.value[k] * w.value[k];
            a -= z.value[2] * (w.gradient[0][0] + w.gradient[1][1]) +
                 w.value[2] * (z.gradient[0][0] + z.gradient[1][1]);
            local(i, j) += a * c.weights[q];
          }
        }
      }
      constraints.distribute_local_to_global(local, b, c.indices, matrix, rhs);
    }
    solve_linear(solution);
    for (const auto i : pressure_dofs)
      solution[i] = 0;
  } else {
    for (unsigned c : {3u, 6u}) {
      auto ids = DoFTools::extract_dofs(dh, fe.component_mask(FEValuesExtractors::Scalar(c)));
      for (auto i : ids)
        solution[i] = 1.;
    }
  }
  auto assemble = [&](double t) {
    matrix = 0;
    rhs = 0;
    for (const auto &c : cache) {
      FullMatrix<double> local(nd, nd);
      Vector<double> b(nd);
      for (unsigned q = 0; q < quadrature.size(); ++q) {
        const auto s = evaluate(c, q, solution), previous = evaluate(c, q, old),
                   f = problem.forcing(c.points[q], t);
        const PointJacobian point(s, previous);
        for (unsigned i = 0; i < nd; ++i) {
          const Field w = c.shapes[q][i];
          b[i] -= residual(p, s, previous, w, f) * c.weights[q];
          for (unsigned j = 0; j < nd; ++j)
            local(i, j) += point.entry(p, c.shapes[q][i], c.shapes[q][j]) * c.weights[q];
        }
      }
      constraints.distribute_local_to_global(local, b, c.indices, matrix, rhs);
    }
    for (const auto &c : faces) {
      FullMatrix<double> local(nd, nd);
      Vector<double> b(nd);
      for (unsigned q = 0; q < c.weights.size(); ++q) {
        Field state, previous;
        for (unsigned i = 0; i < nd; ++i)
          for (unsigned k = 0; k < 7; ++k) {
            state.value[k] += solution[c.indices[i]] * static_cast<Field>(c.shapes[q][i]).value[k];
            previous.value[k] += old[c.indices[i]] * static_cast<Field>(c.shapes[q][i]).value[k];
          }
        const double vn = previous.value[0] * c.normals[q][0] + previous.value[1] * c.normals[q][1];
        for (unsigned i = 0; i < nd; ++i) {
          const Field w = c.shapes[q][i];
          double r = 0;
          for (unsigned k = 0; k < 7; ++k)
            if (k != 2)
              r += .5 * vn * (k < 2 ? p.rho : 1) * state.value[k] * w.value[k];
          // Volume stress is mu*F*F^T; physical outlet traction uses mu*(F*F^T-I).
          if (!p.newtonian())
            r -= p.mu * (c.normals[q][0] * w.value[0] + c.normals[q][1] * w.value[1]);
          b[i] -= r * c.weights[q];
          for (unsigned j = 0; j < nd; ++j)
            for (unsigned k = 0; k < 7; ++k)
              if (k != 2)
                local(i, j) += .5 * vn * (k < 2 ? p.rho : 1) *
                               static_cast<Field>(c.shapes[q][j]).value[k] * w.value[k] *
                               c.weights[q];
        }
      }
      constraints.distribute_local_to_global(local, b, c.indices, matrix, rhs);
    }
  };
  std::ofstream history(out + "/history.csv");
  history << std::setprecision(16)
          << "step,time,newton_iterations,increment_inf,residual_inf,L2_velocity,L2_pressure,L2_F,"
             "energy,min_det_F,log_energy,divergence_L2,weak_divergence_inf,energy_balance_"
             "residual,elapsed_seconds\n";
  std::array<double, 3> spacetime{};
  auto diagnostics = [&](unsigned step, double t, unsigned it, double update,
                         double residual_norm) {
    double mean = 0, volume = 0;
    for (const auto &c : cache)
      for (unsigned q = 0; q < quadrature.size(); ++q) {
        mean += evaluate(c, q, solution).value[2] * c.weights[q];
        volume += c.weights[q];
      }
    mean /= volume;
    double ev = 0, ep = 0, ef = 0, energy = 0, mindet = 1e100, logenergy = 0, divergence = 0,
           balance = 0;
    Vector<double> weak(dh.n_dofs());
    for (const auto &c : cache)
      for (unsigned q = 0; q < quadrature.size(); ++q) {
        const auto s = evaluate(c, q, solution), ex = problem.reference(c.points[q], t);
        const double w = c.weights[q];
        for (unsigned k = 0; k < 2; ++k) {
          ev += std::pow(s.value[k] - ex.value[k], 2) * w;
          energy += .5 * p.rho * s.value[k] * s.value[k] * w;
        }
        ep += std::pow(s.value[2] - mean - ex.value[2], 2) * w;
        for (unsigned k = 3; k < 7; ++k) {
          ef += std::pow(s.value[k] - ex.value[k], 2) * w;
          energy += .5 * p.mu * s.value[k] * s.value[k] * w;
        }
        const double det = determinant(deformation(s));
        mindet = std::min(mindet, det);
        logenergy += det > 0 ? -p.mu * std::log(det) * w : std::numeric_limits<double>::quiet_NaN();
        if (step > 0 && !p.contraction) {
          const auto prev = evaluate(c, q, old), f = problem.forcing(c.points[q], t);
          double b = 0;
          for (unsigned k = 0; k < 2; ++k)
            b += p.rho *
                     (s.value[k] * s.value[k] - prev.value[k] * prev.value[k] +
                      std::pow(s.value[k] - prev.value[k], 2)) /
                     (2 * p.dt) +
                 p.nu * s.gradient[k].norm_square() - f.value[k] * s.value[k];
          for (unsigned k = 3; k < 7; ++k)
            b += p.mu *
                     (s.value[k] * s.value[k] - prev.value[k] * prev.value[k] +
                      std::pow(s.value[k] - prev.value[k], 2)) /
                     (2 * p.dt) +
                 p.mu * p.diffusion_coefficient() * s.gradient[k].norm_square() -
                 p.mu * f.value[k] * s.value[k];
          const auto F = deformation(s);
          b += p.mu * p.mu / (2 * p.lambda) * (Giesekus::stress(F).norm_square() - F.norm_square());
          balance += b * w;
        }
        const double div = s.gradient[0][0] + s.gradient[1][1];
        divergence += div * div * w;
        for (unsigned i = 0; i < nd; ++i)
          weak[c.indices[i]] += div * static_cast<Field>(c.shapes[q][i]).value[2] * w;
      }
    if (p.contraction)
      ev = ep = ef = balance = std::numeric_limits<double>::quiet_NaN();
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    history << step << ',' << t << ',' << it << ',' << update << ',' << residual_norm << ','
            << std::sqrt(ev) << ',' << std::sqrt(ep) << ',' << std::sqrt(ef) << ',' << energy << ','
            << mindet << ',' << logenergy << ',' << std::sqrt(divergence) << ','
            << weak.linfty_norm() << ',' << balance << ',' << elapsed << '\n';
    history.flush();
    if (step > 0) {
      spacetime[0] += p.dt * ev;
      spacetime[1] += p.dt * ep;
      spacetime[2] += p.dt * ef;
    }
    std::cout << "step=" << step << " t=" << t << " Newton=" << it << " |du|inf=" << update
              << " |R|inf=" << residual_norm << " L2(v,p,F)=" << std::sqrt(ev) << ","
              << std::sqrt(ep) << "," << std::sqrt(ef) << " det_min=" << mindet
              << " seconds=" << elapsed << std::endl;
    return mean;
  };
  std::cout << "cells=" << tria.n_active_cells() << " dofs=" << dh.n_dofs()
            << " quadrature_points=" << quadrature.size() << std::endl;
  diagnostics(0, 0, 0, 0, 0);
  double mean = 0;
  for (unsigned step = 1; step <= steps; ++step) {
    old = solution;
    solution_constraints.distribute(solution);
    double update = 0;
    unsigned it = 0;
    bool converged = false;
    for (; it < p.max_newton; ++it) {
      assemble(step * p.dt);
      solve_linear(increment);
      update = increment.linfty_norm();
      solution += increment;
      solution_constraints.distribute(solution);
      if (update < p.tolerance) {
        converged = true;
        ++it;
        break;
      }
    }
    AssertThrow(converged, ExcMessage("Newton failed; run is incomplete"));
    assemble(step * p.dt);
    const double residual_norm = rhs.linfty_norm();
    AssertThrow(std::isfinite(residual_norm) && residual_norm < 1e-9,
                ExcMessage("Post-Newton residual exceeds tolerance"));
    mean = diagnostics(step, step * p.dt, it, update, residual_norm);
  }
  Vector<double> output = solution;
  if (!p.contraction)
    for (const auto i : pressure_dofs)
      output[i] -= mean;
  DataOut<2> data;
  data.attach_dof_handler(dh);
  data.add_data_vector(output, std::vector<std::string>{"velocity_x", "velocity_y", "pressure",
                                                        "F11", "F12", "F21", "F22"});
  data.build_patches(mapping, 2);
  std::ofstream vtk(out + "/solution.vtu");
  data.write_vtu(vtk);
  if (p.contraction) {
    std::ofstream samples(out + "/samples.csv");
    samples << std::setprecision(16) << "x,y,vx,vy,p,F11,F12,F21,F22,stress,detF\n";
    for (const auto &c : cache)
      for (unsigned q = 0; q < c.points.size(); ++q) {
        if (std::abs(c.points[q][0]) > 2 * p.length)
          continue;
        const auto state = evaluate(c, q, solution);
        auto stress = Giesekus::stress(deformation(state));
        stress[0][0] -= 1;
        stress[1][1] -= 1;
        samples << c.points[q][0] << ',' << c.points[q][1];
        for (double value : state.value)
          samples << ',' << value;
        samples << ',' << p.mu * stress.norm() << ',' << determinant(deformation(state)) << '\n';
      }
  }
  auto number = [](double x) {
    std::ostringstream o;
    o << std::setprecision(16) << x;
    return o.str();
  };
  std::ofstream meta(out + "/summary.json");
  meta << std::setprecision(16)
       << "{\n  \"completed\": true,\n  \"cells\": " << tria.n_active_cells()
       << ",\n  \"dofs\": " << dh.n_dofs() << ",\n  \"steps\": " << steps
       << ",\n  \"L2_spacetime_velocity\": "
       << (p.contraction ? "null" : number(std::sqrt(spacetime[0])))
       << ",\n  \"L2_spacetime_pressure\": "
       << (p.contraction ? "null" : number(std::sqrt(spacetime[1])))
       << ",\n  \"L2_spacetime_F\": " << (p.contraction ? "null" : number(std::sqrt(spacetime[2])))
       << ",\n  \"elapsed_seconds\": "
       << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
       << "\n}\n";
}
} // namespace research
