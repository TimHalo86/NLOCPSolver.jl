include("types.jl")


function defineOCP(;
    numStates::Int = 0,
    numControls::Int = 0,
    X0 = fill(NaN,numStates),
    XF = fill(NaN,numStates),
    XL = fill(NaN,numStates),
    XU = fill(NaN,numStates),
    CL = fill(NaN,numControls),
    CU = fill(NaN,numControls))::OCP
    ocp = OCP()
    if numStates < 0 error("States number must be larger than 0") end
    if numControls < 0 error("Control number must be larger than 0") end
    if length(X0) != numStates error("Number of states do not match X0") end
    if length(XF) != numStates error("Number of states do not match XF") end
    if length(XL) != numStates error("Number of states do not match XL") end
    if length(XU) != numStates error("Number of states do not match XU") end
    if length(CL) != numControls error("Number of controls do not match CL") end
    if length(CU) != numControls error("Number of controls do not match CU") end

    ocp.s.states.num    = numStates
    ocp.s.control.num   = numControls
    ocp.b.X0            = X0
    ocp.b.X0_tol        = fill(NaN, numStates)
    ocp.b.XF_tol        = fill(NaN, numStates)
    ocp.b.XL            = XL
    ocp.b.XF            = XF
    ocp.b.XU            = XU
    ocp.b.CL            = CL
    ocp.b.CU            = CU


    return ocp
end

function defineStates!(ocp::OCP, states::Vector{Symbol})
    if length(states) != ocp.s.states.num error("Wrong number of state variable") end
    ocp.s.states.name = states
    return nothing
end

function defineControls!(ocp::OCP, controls::Vector{Symbol})
    if length(controls) != ocp.s.control.num error("Wrong number of state variable") end
    ocp.s.control.name = controls
    return nothing
end
function defineTolerance!(ocp::OCP; X0_tol=fill(NaN, ocp.s.states.num), XF_tol=fill(NaN, ocp.s.states.num))
    if sum(isnan(X0_tol)) != ocp.s.states.num
        ocp.s.X0slack == true
        ocp.b.X0_tol =X0_tol
    end
    if sum(isnan(XF_tol)) != ocp.s.states.num
        ocp.s.XFslack == true
        ocp.b.XF_tol =XF_tol
    end
    return nothing
end


function ValidateScheme(scheme)
    num = maximum(size(scheme))
    for i in 1:num
        if scheme[i] ∉ [:bkwEuler, :trapezoidal, :RK1, :RK2, :RK3, :RK4]
            return false, i
        end
    end
    return true
end


function ConfigurePredefined(ocp::OCP; kwargs...)::OCPFormulation
    # all inputs here: tfDV, tf, Np, Integration Scheme
    # Integration scheme: bkwEuler, trapezoidal
    # Trajectory methods: MultipleShooting Collocation
    OCPForm = OCPFormulation()
    OCPForm.mdl = JuMP.Model()
    set_silent(OCPForm.mdl)
    kw = Dict(kwargs)
    if haskey(kw, :tfDV)
        OCPForm.tfDV = get(kw, :tfDV, 0)
    end
    if haskey(kw, :tf)
        if OCPForm.tfDV == true
            @warn "Conflict intructions to put tf as design variable, change to fixed horizon"
            OCPForm.tfDV = false
        end
        OCPForm.tf = get(kw, :tf, 0)
        if OCPForm.tf <= ocp.b.tfMin || OCPForm.tf >= ocp.b.tfMax
            error("Please make sure tf ∈ [$(ocp.b.tfMin), $(ocp.b.tfMax)]")
        end
    else
        if OCPForm.tfDV == false
            error("Either treat tf as design variable or use fixed time horizon")
        end
    end

    if haskey(kw, :Np)
        Np = get(kw, :Np, 0)
        NpType = typeof(Np)
        if NpType != Int64 
            if NpType <: Real
                @warn "Round Np to nearest integer"
                OCPForm.Np = Int(round(a))
            else
                error("Wrong type of number of points, should be Int64")
            end
        else
            OCPForm.Np = Np
        end
        if OCPForm.Np <= 0
            error("Make sure Np is an integer and is larger than 0")
        end
    else
        error("No number of point input")
    end

    if haskey(kw, :IntegrationScheme)
        IntegrationScheme_single = get(kw, :IntegrationScheme, 0)
        if IntegrationScheme_single ∈ [:bkwEuler, :trapezoidal, :RK1, :RK2, :RK3, :RK4]
            OCPForm.IntegrationScheme = vec(repeat([IntegrationScheme_single],OCPForm.Np-1))
        else
            @warn "$IntegrationScheme is not implemented, use default (bkwEuler)"
        end
    end

    if OCPForm.tfDV == false
        OCPForm.tw = 1 / (OCPForm.Np - 1) .* ones(OCPForm.Np - 1)
        OCPForm.TInt = OCPForm.tf .*  OCPForm.tw
    elseif OCPForm.tfDV == true
        OCPForm.tw = 1 / (OCPForm.Np - 1) .* ones(OCPForm.Np - 1)
        OCPForm.tf = @variable(OCPForm.mdl, ocp.b.tfMin <= tf <= ocp.b.tfMax)
        OCPForm.TInt = @expression(OCPForm.mdl, [idx = 1:OCPForm.Np - 1], OCPForm.tf * OCPForm.tw[idx])
    end

    if !haskey(kw, :dx)
        error("No dynamics here")
    else
        OCPForm.dx = Vector{Any}(nothing, OCPForm.Np)
        OCPForm.dx[1:OCPForm.Np] .= get(kw, :dx, 0)
    end
    if haskey(kw, :params)
        params = get(kw, :params, 0)
        if size(params, 1) == OCPForm.Np
            OCPForm.params = Matrix{Any}(undef, size(params, 1), size(params, 2))
            OCPForm.params[1:size(params, 1), 1:size(params, 2)] = params
        else
            if size(params, 2) == 1 && size(params, 1) != 1
                OCPForm.params = repeat(params', OCPForm.Np)
            elseif size(params, 1) == 1 && size(params, 2) != 1
                OCPForm.params = repeat(params, OCPForm.Np)
            elseif size(OCPForm.params, 1) == 1 && size(OCPForm.params, 1) == 1
                OCPForm.params = repeat(params, OCPForm.Np)
            else
                error("Size of param ($(size(params, 1)) * $(size(params, 2)) is not expected)")
            end
        end
    end
    if size(OCPForm.params, 1) == 0 || size(OCPForm.params, 2) == 0
        OCPForm.params = Matrix{Any}(undef, OCPForm.Np, 1)
    end

    OCPForm.cons = Vector{Any}(nothing, OCPForm.Np)
    if haskey(kw, :cons)
        OCPForm.cons .= get(kw, :cons, 0)
    end

    OCPForm.expr = Vector{Any}(nothing, OCPForm.Np)

    if haskey(kw, :expr)
        OCPForm.expr .= get(kw, :expr, 0)
    end
    return OCPForm
end

function CheckOCPFormulation(ocp::OCP, OCPForm::OCPFormulation)
    tw = OCPForm.tw
    if abs(1 - sum(tw)) > 1e-4
        error("Wrong weights of tf")
    end

    if length(OCPForm.tw) != length(OCPForm.TInt) error("Size of tw and TInt do not match") end
    if (OCPForm.Np < 0) || (typeof(OCPForm.Np) != Int64) error("Wrong input of Np") end
    if OCPForm.tfDV == true && typeof(OCPForm.tf) != JuMP.VariableRef error("Wrong type of tf since tfDV is $(OCPForm.tfDV)") end
    if OCPForm.tfDV == false 
        if !(typeof(OCPForm.tf) <: Real)
            error("Wrong type of tf since tfDV is $(OCPForm.tfDV)") 
        else
            if OCPForm.tf < ocp.b.tfMin || OCPForm.tf > ocp.b.tfMax
                error("Please make sure tf ∈ [$(ocp.b.tfMin), $(ocp.b.tfMax)]")
            end
        end
    end

    for i in eachindex(OCPForm.TInt)
        if OCPForm.tfDV == true 
            if OCPForm.TInt[i] != OCPForm.tw[i] * OCPForm.tf
                error("Weight and TInt do not correspond")
            end
        else
            if abs(OCPForm.TInt[i] - OCPForm.tw[i] * OCPForm.tf) > 1e-4
                error("Weight and TInt do not correspond")
            end
        end
    end
    ValidationInfo = ValidateScheme(OCPForm.IntegrationScheme)
    if length(ValidationInfo) == 2 # Midpoint Collocation?
        error("$(OCPForm.IntegrationScheme[ValidationInfo[2]]) is not implemented")
    end

    if size(OCPForm.params, 1) != OCPForm.Np
        error("Size [$(size(OCPForm.params, 1)), $(size(OCPForm.params, 2))] do not match number of points $(OCPForm.Np) ")
    end

    return nothing
end


function defineSolver!(OCPForm::OCPFormulation, SolverName::Symbol, Options::Tuple)
    if solver_name(OCPForm.mdl) == "No optimizer attached."
        if SolverName == :Ipopt
            set_optimizer(OCPForm.mdl, Ipopt.Optimizer)
            set_attributes(OCPForm.mdl, Options...)
        else
            error("Solver $SolverName is not implemented")
        end
    else
        @warn "$(solver_name(OCPForm.mdl)) is set, can not change to $SolverName"
    end
    return nothing
end


function CalXvar(Scheme)
    count = 0
    for i in 1:maximum(size(Scheme))
        if Scheme[i] == :RK2
            count += 1
        elseif Scheme[i] == :RK3
            count += 2
        elseif Scheme[i] == :RK4
            count += 3
        end
    end
    return Int32(ceil(count))
end

function OCPdef!(ocp::OCP, OCPForm::OCPFormulation)
    ocp.f = OCPForm
    CheckOCPFormulation(ocp, OCPForm)
    ## TODO Configure the model setting
    if solver_name(OCPForm.mdl) == "No optimizer attached."
        defineSolver!(OCPForm, ocp.s.solver.name, ocp.s.solver.settings)
    end
    ## TODO Write the Lower and upper bound
    ocp.s.states.pts = ocp.s.control.pts = OCPForm.Np
    ocp.p.x = @variable(OCPForm.mdl, ocp.b.XL[i] <= x[j in 1:ocp.s.states.pts, i in 1:ocp.s.states.num] <= ocp.b.XU[i])
    ocp.p.u = @variable(OCPForm.mdl, ocp.b.CL[i] <= u[j in 1:ocp.s.control.pts, i in 1:ocp.s.control.num] <= ocp.b.CU[i])
    OCPForm.mdl[:x] = ocp.p.x
    OCPForm.mdl[:u] = ocp.p.u

    ## fix the initial states
    
    if !ocp.s.X0slack
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.X0[st]) fix(OCPForm.mdl[:x][1, st], ocp.b.X0[st]; force = true) end
        end
    else
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.X0[st])
                if !isnan(ocp.b.X0_tol[st])
                    @constraint(OCPForm.mdl, ocp.b.X0[st] - ocp.b.X0_tol[st] <= OCPForm.mdl[:x][1, st] <= ocp.b.X0[st] + ocp.b.X0_tol[st])
                else
                    @constraint(OCPForm.mdl, OCPForm.mdl[:x][1, st] == ocp.b.X0[st])
                end
            end
        end
    end
    
    if !ocp.s.XFslack
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.XF[st]) fix(OCPForm.mdl[:x][end, st], ocp.b.XF[st]; force = true) end
        end
    else
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.XF[st])
                if !isnan(ocp.b.XF_tol[st])
                    @constraint(OCPForm.mdl, ocp.b.XF[st] - ocp.b.XF_tol[st] <= OCPForm.mdl[:x][end, st] <= ocp.b.XF[st] + ocp.b.XF_tol[st])
                else
                    @constraint(OCPForm.mdl, OCPForm.mdl[:x][end, st] == ocp.b.XF[st])
                end
            end
        end
    end
    
    # Parameters
    Nonparam = Any
    if isassigned(OCPForm.params, 1)
        ocp.p.params = @variable(OCPForm.mdl, Params[i = 1:size(OCPForm.params, 1), j = 1:size(OCPForm.params, 2)] in Parameter(OCPForm.params[i, j]) )
    end





    # Dynamical Constraints
    δx = Matrix{Any}(undef, ocp.s.states.pts, ocp.s.states.num)

    ColumnCounter = 1

    numXvar = CalXvar(OCPForm.IntegrationScheme)
    if numXvar > 0
        ocp.p.xvar = @variable(OCPForm.mdl, xvar[j in 1:numXvar, i in 1:ocp.s.states.num])
    end



    for j in 1:ocp.s.states.pts - 1
        if isassigned(OCPForm.params, j)
            param = ocp.p.params[j, :]
        else
            param = Nonparam
        end
        if OCPForm.IntegrationScheme[j] == :RK2
            xk1 = xvar[ColumnCounter, :]
            ColumnCounter += 1
            k1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            k2 = @expression(OCPForm.mdl, OCPForm.dx[j](xk1, ocp.p.u[j, :], param))
            @constraint(OCPForm.mdl, [i = 1:ocp.s.states.num],  xk1[i] - ocp.p.x[j, i] == k1[i] * OCPForm.TInt[j])
            δx[j, :] = @expression(OCPForm.mdl, k1 ./ 2 .+ k2 ./ 2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])
        elseif OCPForm.IntegrationScheme[j] == :RK3
            xk1 = xvar[ColumnCounter,:]
            xk2 = xvar[ColumnCounter+1,:]
            ColumnCounter += 2

            k1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            k2 = @expression(OCPForm.mdl, OCPForm.dx[j](xk1, ocp.p.u[j, :], param))
            k3 = @expression(OCPForm.mdl, OCPForm.dx[j](xk2, ocp.p.u[j, :], param))
            
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk1[i] - ocp.p.x[j, i] == k1[i] * OCPForm.TInt[j]/2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk2[i] - ocp.p.x[j, i] == -k1[i] * OCPForm.TInt[j] + 2*k2[i] * OCPForm.TInt[j])

            δx[j, :] = @expression(OCPForm.mdl, (k1 .+  4*k2 .+ k3) ./ 6 )
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        elseif OCPForm.IntegrationScheme[j] == :RK4
            xk1 = xvar[ColumnCounter,:]
            xk2 = xvar[ColumnCounter+1,:]
            xk3 = xvar[ColumnCounter+2,:]
            ColumnCounter += 3

            k1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            k2 = @expression(OCPForm.mdl, OCPForm.dx[j](xk1, ocp.p.u[j, :], param))
            k3 = @expression(OCPForm.mdl, OCPForm.dx[j](xk2, ocp.p.u[j, :], param))
            k4 = @expression(OCPForm.mdl, OCPForm.dx[j](xk3, ocp.p.u[j, :], param))
            
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk1[i] - ocp.p.x[j, i] == k1[i] * OCPForm.TInt[j]/2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk2[i] - ocp.p.x[j, i] == k2[i] * OCPForm.TInt[j]/2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk3[i] - ocp.p.x[j, i] == k3[i] * OCPForm.TInt[j])

            δx[j, :] = @expression(OCPForm.mdl, (k1 .+  2*k2 .+ 2*k3 .+ k4) ./ 6 )
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        elseif OCPForm.IntegrationScheme[j] == :RK1
            δx[j, :] = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        elseif OCPForm.IntegrationScheme[j] == :bkwEuler
            δx[j, :] = @expression(OCPForm.mdl, OCPForm.dx[j + 1](ocp.p.x[j + 1, :], ocp.p.u[j + 1, :], param))
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] == ocp.p.x[j, i] + δx[j, i] * OCPForm.TInt[j] )

        elseif OCPForm.IntegrationScheme[j] == :trapezoidal
            δx1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            δx2 = @expression(OCPForm.mdl, OCPForm.dx[j + 1](ocp.p.x[j + 1, :], ocp.p.u[j + 1, :], param))
            δx[j, :] = @expression(OCPForm.mdl, δx1 ./ 2 .+ δx2 ./ 2 )
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        end

    end

    ocp.p.δx = δx
  
    
    # Inner States Constraints
    for j in 1:ocp.s.states.pts 
        if isassigned(OCPForm.params, j)
            param = ocp.p.params[j, :]
        else
            param = Nonparam
        end

        if !isnothing(OCPForm.cons[j])
            @constraints(OCPForm.mdl, begin OCPForm.cons[j](OCPForm.mdl[:x][j, :], OCPForm.mdl[:u][j, :], param) >= 0 end)
        end
    end

    ocp.p.tV = [0.0; cumsum(OCPForm.TInt)]
    ocp.f = OCPForm
    return nothing

end

