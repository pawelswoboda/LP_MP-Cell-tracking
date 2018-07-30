
#include "cell_tracking.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

int main(int argc, char** argv) {
Solver<LP<FMC_CELL_TRACKING_FLOW>,StandardVisitor> solver(argc,argv);
auto input = cell_tracking_parser_2d::parse_file(solver.get_input_file());
auto& c = solver.template GetProblemConstructor<0>();
c.construct(input);
return solver.Solve();
}