#!/usr/bin/python

from collections import namedtuple
solver = namedtuple("solver_info", "preamble FMC LP parse_fun filename")

preamble = """
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
"""

solvers = [
    solver(preamble, 'FMC_CELL_TRACKING_MOTHER_MACHINE', 'LP', 'cell_tracking_parser_mother_machine::ParseProblemMotherMachine', "cell_tracking_mother_machine.cpp"),
    solver(preamble, 'FMC_CELL_TRACKING', 'LP', 'cell_tracking_parser_2d::ParseProblem', "cell_tracking.cpp"),
    solver(preamble, 'FMC_CELL_TRACKING_WITH_DIVISION_DISTANCE', 'LP', 'cell_tracking_parser_2d::parse_problem_with_division_distance', "cell_tracking_with_division_distance.cpp"),
    solver(preamble, 'FMC_CELL_TRACKING_DUPLICATE_EDGES', 'LP', 'cell_tracking_parser_2d::ParseProblem', "cell_tracking_duplicate_edges.cpp"),
    solver(preamble, 'FMC_CELL_TRACKING_DIVISION_DISTANCE_DUPLICATE_EDGES', 'LP', 'cell_tracking_parser_2d::parse_problem_with_division_distance', "cell_tracking_division_distance_duplicate_edges.cpp"),
    solver(preamble, 'FMC_CELL_TRACKING_FINE_DECOMPOSITION', 'LP', 'cell_tracking_parser_2d::ParseProblem', "cell_tracking_fine_decomposition.cpp")
    ]

for e in solvers:
   f = open(e.filename,'w')
   f.write(e.preamble)
   f.write("\nusing namespace LP_MP;\n\nint main(int argc, char** argv) {\nSolver<")
   f.write(e.FMC + "," + e.LP + ",StandardVisitor> solver(argc,argv);\n")
   f.write("solver.ReadProblem(" + e.parse_fun + "<Solver<" + e.FMC + "," + e.LP + ",StandardVisitor>>);\n")
   f.write("return solver.Solve();\n}")
   f.close()
