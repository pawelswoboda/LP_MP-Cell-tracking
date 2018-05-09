
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

int main(int argc, char** argv) {
Solver<FMC_CELL_TRACKING_DIVISION_DISTANCE_DUPLICATE_EDGES,LP<FMC_CELL_TRACKING_DIVISION_DISTANCE_DUPLICATE_EDGES>,StandardVisitor> solver(argc,argv);
solver.ReadProblem(cell_tracking_parser_2d::parse_problem_with_division_distance<Solver<FMC_CELL_TRACKING_DIVISION_DISTANCE_DUPLICATE_EDGES,LP<FMC_CELL_TRACKING_DIVISION_DISTANCE_DUPLICATE_EDGES>,StandardVisitor>>);
return solver.Solve();
}