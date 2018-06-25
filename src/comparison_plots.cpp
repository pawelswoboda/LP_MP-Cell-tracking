#include "cell_tracking.h"
#include "solver.hxx"
#include "visitors/sqlite_visitor.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
    using solver_type = Solver<FMC_CELL_TRACKING,LP<FMC_CELL_TRACKING>,SqliteVisitor<StandardVisitor>>;

    std::vector<std::string> options = 
    {
        {""},
        {"-i"}, {"/BS/discrete_opt/work/datasets/cell-tracking/20170607_Berkeley_FullDataset_GarciaLab_MULTI.5.7.85/tr2d_problem.jug"},
        {"--primalComputationStart"}, {"10000"},
        {"--maxIter"}, {"2000"},
        {"--lowerBoundComputationInterval"}, {"10"},
        {"--databaseFile"},{"cell_tracking_plots.db"},
        {"--datasetName"},{"Garcia"},
        {"--algorithmFMC"},{"CELL_TRACKING"}
    };

    { 
        auto srmp_options = options;
        srmp_options.push_back("--algorithmName");
        srmp_options.push_back("srmp");
        solver_type solver(srmp_options);
        solver.ReadProblem(cell_tracking_parser_2d::ParseProblem<solver_type>);
        solver.Solve();
    }

    {
        auto damped_uniform_options = options;
        damped_uniform_options[6] = "100";
        damped_uniform_options.push_back("--standardReparametrization"); 
        damped_uniform_options.push_back("damped_uniform"); 
        damped_uniform_options.push_back("--algorithmName");
        damped_uniform_options.push_back("mcf");
        Solver<FMC_CELL_TRACKING_FLOW,LP<FMC_CELL_TRACKING_FLOW>, SqliteVisitor<StandardVisitor>> solver(damped_uniform_options);
        solver.Solve();
    }

    exit(0);

    { 
        auto partition_options = options;
        partition_options.push_back("--reparametrizationType");
        partition_options.push_back("partition");
        partition_options.push_back("--innerIteration");
        partition_options.push_back("5");
        partition_options.push_back("--algorithmName");
        partition_options.push_back("partition,5");
        solver_type solver(partition_options);
        solver.ReadProblem(cell_tracking_parser_2d::ParseProblem<solver_type>);
        solver.Solve();
    }

    { 
        auto partition_options = options;
        partition_options.push_back("--reparametrizationType");
        partition_options.push_back("partition");
        partition_options.push_back("--innerIteration");
        partition_options.push_back("10");
        partition_options.push_back("--algorithmName");
        partition_options.push_back("partition,10");
        solver_type solver(partition_options);
        solver.ReadProblem(cell_tracking_parser_2d::ParseProblem<solver_type>);
        solver.Solve();
    }

    { 
        auto partition_options = options;
        partition_options.push_back("--reparametrizationType");
        partition_options.push_back("overlapping_partition");
        partition_options.push_back("--innerIteration");
        partition_options.push_back("5");
        partition_options.push_back("--algorithmName");
        partition_options.push_back("overlapping partition,5");
        solver_type solver(partition_options);
        solver.ReadProblem(cell_tracking_parser_2d::ParseProblem<solver_type>);
        solver.Solve(); 
    }

    { 
        auto partition_options = options;
        partition_options.push_back("--reparametrizationType");
        partition_options.push_back("overlapping_partition");
        partition_options.push_back("--innerIteration");
        partition_options.push_back("10");
        partition_options.push_back("--algorithmName");
        partition_options.push_back("overlapping partition,10");
        solver_type solver(partition_options);
        solver.ReadProblem(cell_tracking_parser_2d::ParseProblem<solver_type>);
        solver.Solve(); 
    }

}
