# Author: Siyuan Yu and Congkai Shen 
# Position: Graduate Student in University of Michigan
# Date Created: 02/28/2024


using JuMP
using Parameters
using DataFrames
using Ipopt

export  States,
        Control,
        Solver,
        OCPFormulation,
        OCPSetting,
        OCPBound,
        OCPParameter,
        OCPResults,
        OCP


_Ipopt_MPC_defaults = (
    "mu_strategy" => "adaptive",
    "max_iter" => 100,
    "tol" => 6e-3,
    "dual_inf_tol" => 2.,
    "constr_viol_tol" => 3e-1,
    "compl_inf_tol" => 3e-1,
    "acceptable_tol" => 2e-2,
    "acceptable_constr_viol_tol" => 0.02,
    "acceptable_dual_inf_tol" => 1e10,
    "acceptable_compl_inf_tol" => 0.02,
    "warm_start_init_point" => "yes",
    "fixed_variable_treatment" => "relax_bounds",
    "max_cpu_time" => 0.1,
    "print_level" => 0
    )

_Ipopt_defaults = (
    "mu_strategy" => "adaptive",
    "max_iter" => 1000,
    "warm_start_init_point" => "no",
    "fixed_variable_treatment" => "relax_bounds",
    "max_cpu_time" => 10.0,
    "print_level" => 0
    )

@with_kw mutable struct States
    num::Int                                    = 0                                 # Number of States
    name::Vector{Symbol}                        = Vector{Symbol}[]                  # Names of States Variable
    pts::Int                                    = 0                                 # Number of points in States Variable
end

@with_kw mutable struct Control
    num::Int                                    = 0                                 # Number of Control
    name::Vector{Symbol}                        = Vector{Symbol}[]                  # Names of Constrol Variable
    pts::Int                                    = 0                                 # Number of points in Control Variable
end

@with_kw mutable struct Solver
    name::Symbol                                = :Ipopt                            # Default solver name
    settings::Tuple                             = _Ipopt_defaults                   # Default settings
end

@with_kw mutable struct OCPFormulation{ T <: Number }
    tfDV::Bool                                  = false                             # Determines whether tf is a design variable
    Np::Int64                                   = 0                                 # NckPoint 
    IntegrationScheme::Vector{Symbol}           = Vector{Symbol}()                         # Integration integrationScheme
    tf::Any                                     = Any                               # Total time, type depend on the 
    tw::Vector{Float64}                         = Vector{Float64}()                 # Time weight
    TInt::Vector{Any}                           = Vector{Any}[]                     # Depends on terminal constraints (Nonlinear Expr) or fixed time horizon (Float64)
    mdl::JuMP.Model                             = JuMP.Model()                      # JuMP model
    dx::Vector{Any}                             = Vector{Any}()                     # Dynamics function handle here
    cons::Vector{Any}                           = Vector{Any}()                     # Constraints function handle here
    expr::Vector{Any}                           = Vector{Any}()                     # Expression function handle here
    params::Matrix{Any}                         = Matrix{Any}(undef, 0, 0)                     # Place Holder in case user uses 
end

@with_kw mutable struct OCPSetting{ T <: Number }
    states::States                              = States()                          # States structure in OCPSetting
    control::Control                            = Control()                         # Control structure in OCPSetting
    solver::Solver                              = Solver()                          # Solver structure in OCPSetting
    InternalLogging::Bool                       = true                              # Bool for logging data internally
    X0slack::Bool                               = false                             # Boolean for using tolerance on initial states
    XFslack::Bool                               = false                             # Boolean for using tolerance on final states
end

@with_kw mutable struct OCPBound{ T <: Number }
    tfMin::Float64                              = 0.0                               # tf minimum
    tfMax::Float64                              = 100                               # tf maximum
    X0::Vector{T}                               = Vector{T}[]                       # Initial Condition
    X0_tol::Vector{T}                           = Vector{T}[]                       # Initial Condition Slack Variable 
    XF::Vector{T}                               = Vector{T}[]                       # Final Condition Slack Variable 
    XF_tol::Vector{T}                           = Vector{T}[]                       # Final Condition Slack Variable 
    XL::Vector{T}                               = Vector{T}[]                       # Lower bound of states variable
    XU::Vector{T}                               = Vector{T}[]                       # Upper bound of states variable
    CL::Vector{T}                               = Vector{T}[]                       # Lower bound of control variable
    CU::Vector{T}                               = Vector{T}[]                       # Upper bound of control variable
end

@with_kw mutable struct OCPParameter{ T <: Number }
    x::Matrix{Any}                              = Matrix{Any}(undef,0,0)            # Holder for JuMP nonlinear variable(collocation); nonlinear expression (single shooting)
    u::Matrix{VariableRef}                      = Matrix{VariableRef}(undef,0,0)    # Control inputs are always variable references
    params::Matrix{VariableRef}                 = Matrix{VariableRef}(undef, 0, 0)  # paramters used in dynamics function
    tV::Any                                     = Any                               # Time point
    xvar::Matrix{VariableRef}                   = Matrix{VariableRef}(undef,0, 0)   # Used for register variable states in multiple shooting method (Not used in collocation method)
    δx::Matrix{Any}                             = Matrix{Any}(undef,0,0)            # Place Holder for derivatives
end

@with_kw mutable struct OCPResults{ T <: Number }
    X::Matrix{Float64}                          = Matrix{Float64}(undef,0,0)        # State variable value
    U::Matrix{Float64}                          = Matrix{Float64}(undef,0,0)        # Control variable value
    Tst::Vector{Float64}                        = Vector{Float64}()                 # Time value
    dt::Vector{Float64}                         = Vector{Float64}()                 # Time Interval Value
    Status::Symbol                              = :InFeasible                       # Status symbol: :Infeasible, :UserLimit, :Optimal
    IterNum::Int64                              = 0                                 # Iteration number from Solver
    EvalNum::Int64                              = 0                                 # Evaluation number from Solver
    TerminalStatus::MOI.TerminationStatusCode   = MOI.OTHER_ERROR                   # Math operation interface termination status
    TSolve::Float64                             = 0.0                               # Solve time
    Objval::Float64                             = 0.0                               # objective value
    Dfs::Vector{DataFrame}                      = Vector{DataFrame}()               # DataFrame results of OCP
end

@with_kw mutable struct OCP{ T <: Number }
    s::OCPSetting{T}                            = OCPSetting{T}()
    b::OCPBound{T}                              = OCPBound{T}()
    f::OCPFormulation{T}                        = OCPFormulation{T}()
    p::OCPParameter{T}                          = OCPParameter{T}()    
    r::OCPResults{T}                            = OCPResults{T}()
end
OCPFormulation() = OCPFormulation{Float64}()
OCP() = OCP{Float64}()
