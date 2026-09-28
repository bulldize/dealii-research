#include <deal.II/base/parameter_handler.h>
#include <iostream>
#include <research/manufactured.h>
#include <research/solver.h>
int main(int argc, char **argv) {
  try {
    dealii::ParameterHandler h;
    h.declare_entry("subdivisions", "16", dealii::Patterns::Integer(2));
    h.declare_entry("dt", "0.001", dealii::Patterns::Double(1e-10));
    h.declare_entry("end_time", "0.1", dealii::Patterns::Double(1e-10));
    h.declare_entry("rho", "1", dealii::Patterns::Double(1e-10));
    h.declare_entry("nu", "1", dealii::Patterns::Double(1e-10));
    h.declare_entry("mu", "1", dealii::Patterns::Double(1e-10));
    h.declare_entry("lambda", "1", dealii::Patterns::Double(0));
    h.declare_entry("newton_tolerance", "1e-12", dealii::Patterns::Double(1e-15));
    h.declare_entry("max_newton", "15", dealii::Patterns::Integer(1));
    h.declare_entry("problem", "manufactured",
                    dealii::Patterns::Selection("manufactured|contraction"));
    h.declare_entry("mesh_file", "", dealii::Patterns::Anything());
    h.declare_entry("diffusion", "-1", dealii::Patterns::Double(-1));
    h.declare_entry("output_interval", "0", dealii::Patterns::Integer(0));
    if (argc != 3) {
      std::cerr << "Usage: giesekus config.prm output_directory\n";
      return 2;
    }
    h.parse_input(argv[1]);
    research::Parameters p;
    p.subdivisions = h.get_integer("subdivisions");
    p.dt = h.get_double("dt");
    p.end = h.get_double("end_time");
    p.rho = h.get_double("rho");
    p.nu = h.get_double("nu");
    p.mu = h.get_double("mu");
    p.lambda = h.get_double("lambda");
    p.tolerance = h.get_double("newton_tolerance");
    p.max_newton = h.get_integer("max_newton");
    p.contraction = h.get("problem") == "contraction";
    p.mesh_file = h.get("mesh_file");
    p.diffusion = h.get_double("diffusion");
    p.output_interval = h.get_integer("output_interval");
    if (!p.contraction && p.newtonian())
      throw std::runtime_error("lambda=0 is only supported for the contraction Newtonian limit");
    research::Problem problem{
        "paper-section-5.2", research::exact,
        [p](const dealii::Point<2> &x, double t) { return research::source(x, t, p); }};
    if (p.contraction) {
      problem.name = "paper-section-5.3";
      problem.reference = [](const dealii::Point<2> &, double) {
        research::Field s;
        s.value[3] = s.value[6] = 1;
        return s;
      };
      problem.forcing = [](const dealii::Point<2> &, double) { return research::Field{}; };
    }
    research::run(p, problem, argv[2]);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
