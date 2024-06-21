include("setup.jl")


function OptSolve!(ocp::OCP)
    JuMP.optimize!(ocp.f.mdl)
    ocp.r.TerminalStatus = termination_status(ocp.f.mdl)
    ocp.r.Status = RetrieveSolveStatus(ocp.r.TerminalStatus)
    ocp.r.Objval = objective_value(ocp.f.mdl)
    ocp.r.TSolve = solve_time(ocp.f.mdl)
    ocp.r.IterNum = barrier_iterations(ocp.f.mdl)
    ocp.r.EvalNum = ocp.r.EvalNum + 1
    GetOptimizeValue!(ocp)
    return nothing
end



function ExprIntegral(ocp::OCP)
    cost = @expression(ocp.f.mdl, 0)
    Nonparam = Any
    for j in 2:ocp.f.Np
        if !isnothing(ocp.f.expr[j])
            if isassigned(ocp.f.params, j)
                param = ocp.p.params[j, :]
            else
                param = Nonparam
            end

            if ocp.f.IntegrationScheme ∈ [:RK1, :RK2, :RK3, :RK4]
                cost = @expression(ocp.f.mdl, cost + ocp.f.expr[j](ocp.p.x[j - 1, :], ocp.p.u[j - 1 - 1, :], param) * ocp.f.TInt[j - 1])
            elseif ocp.f.IntegrationScheme == :trapezoidal
                cost = @expression(ocp.f.mdl, cost + ocp.f.expr[j](ocp.p.x[j, :], (ocp.p.u[j - 1, :] + ocp.p.u[j, :]) / 2, param) * ocp.f.TInt[j - 1])
            elseif ocp.f.IntegrationScheme == :bkwEuler
                cost = @expression(ocp.f.mdl, cost + ocp.f.expr[j](ocp.p.x[j, :], ocp.p.u[j, :], param) * ocp.f.TInt[j - 1])
            end
        end
    end
    return cost
end




function RetrieveSolveStatus(status::MOI.TerminationStatusCode)
    SolvingStatus = [:Optimal, :UserLimit, :InFeasible]
    OptimalList = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
    LimitList = [MOI.ITERATION_LIMIT, MOI.TIME_LIMIT, MOI.NODE_LIMIT, MOI.SOLUTION_LIMIT, MOI.MEMORY_LIMIT, MOI.OBJECTIVE_LIMIT, MOI.NORM_LIMIT, MOI.OTHER_LIMIT ]
    
    if status ∈ OptimalList
        return SolvingStatus[1]
    elseif status ∈ LimitList
        return SolvingStatus[2]
    else
        return SolvingStatus[3]
    end
end

function GetOptimizeValue!(ocp::OCP)
    if ocp.r.Status ∈ [:Optimal, :UserLimit]
        ocp.r.X = value.(ocp.p.x)
        ocp.r.U = value.(ocp.p.u)
        if ocp.f.tfDV
            ocp.r.Tst = value.(ocp.p.tV)
            ocp.r.dt = value.(ocp.f.TInt)
        else
            ocp.r.Tst = ocp.p.tV
            ocp.r.dt = ocp.f.TInt
        end
        dfs = ResultsToDataFrame(ocp)
        if ocp.s.InternalLogging
            push!(ocp.r.Dfs, dfs)
        else
            ocp.r.Dfs = Vector{DataFrame}()
            push!(ocp.r.Dfs, dfs)
        end
    else
        @warn "Solution is infeasible"
    end



end

function ResultsToDataFrame(ocp::OCP)
    dfs = DataFrame()
    dfs[!, :t] = ocp.r.Tst
    for st in 1:ocp.s.states.num
        dfs[!, ocp.s.states.name[st]] = ocp.r.X[:,st]
    end
    for ctr in 1:ocp.s.control.num
        dfs[!, ocp.s.control.name[ctr]] = ocp.r.U[:,ctr]
    end
    return dfs
end


function CreateEmptyFormulation()::OCPFormulation
    return OCPFormulation()
end

function DeleteElement(ConsExpr, index)
    ConsExpr[index] .= nothing
end


function WarmStart(ocp::OCP)
    flag = false
    try get_attribute(ocp.f.mdl, "warm_start_init_point")
      if get_attribute(ocp.f.mdl, "warm_start_init_point") == "yes"
        flag = true
      else
        flag = false
      end
    catch
      flag = false
    end
    if flag == true
      set_start_value.(ocp.p.x, [ocp.b.X0'; ocp.r.X[2:end, :]])
      set_start_value.(ocp.p.u, [ocp.r.U[2:end, :]; ocp.r.U[end, :]'])
    end
    return nothing
end
  

function UpdateX0!(ocp::OCP, X0)
    ocp.b.X0 = X0
    for i = 1:1:ocp.s.states.num
        fix(ocp.p.x[1, i], ocp.b.X0[i]; force = true)
    end
    WarmStart(ocp)
end