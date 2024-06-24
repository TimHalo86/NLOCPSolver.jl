isdefined(Base, :__precompile__) && __precompile__()
module NLOCPSolver

using JuMP
using Parameters
using DataFrames
using Ipopt

include("types.jl")

include("utils.jl")

export  OptSolve!,
        ExprIntegral,
        RetrieveSolveStatus,
        GetOptimizeValue!,
        ResultsToDataFrame,
        CreateEmptyFormulation,
        DeleteElement,
        WarmStart,
        UpdateX0!

include("setup.jl")

export  defineOCP,
        defineStates!,
        defineControls!,
        defineTolerance!,
        ValidateScheme,
        ConfigurePredefined,
        CheckOCPFormulation,
        defineSolver!,
        CalXvar,
        OCPdef!

end