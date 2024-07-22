using JuMP
using NLOCPSolver
using Test

la = 1.56 
lb = 1.64
m = 2020
g = 9.81
Izz = 4095
h = 0.6
mu = 0.8

function MagicFormula(alpha, Fz, mu)
    B =  5.68   # Input Q2b value here
    C =  1.817   # Input Q2b value here
    Fy =  mu*Fz*sin(C*atan(B/mu*alpha))  # Lateral force calculation
    return Fy
end

function bicycleModel_expr(states, controls, parameters)
    x = states[1]
    y = states[2]
    v = states[3]
    r = states[4]
    ψ = states[5]
    ux = states[6]
    δf = states[7]
    ax = controls[1]
    dδf = controls[2]
    Fzf = m*g*la/(la+lb) - m*h/(la+lb)*ax # Front axle load
    Fzr = m*g*lb/(la+lb) + m*h/(la+lb)*ax # Rear axle load

    αf = δf - (v+la*r)/ux # Front slip angle
    αr = -(v-lb*r)/ux # Rear slip angle

    Fyf = MagicFormula(αf, Fzf, mu) # Front lateral force
    Fyr = MagicFormula(αr, Fzr, mu) # Rear lateral force

    dstates = Vector{Any}(undef,7)
    dstates[1]         = ux*cos(ψ) - v*sin(ψ) 
    dstates[2]         = ux*sin(ψ) + v*cos(ψ)  
    dstates[3]         = (Fyf+Fyr)/m - ux*r  
    dstates[4]         = (Fyf*la-Fyr*lb)/Izz 
    dstates[5]         = r
    dstates[6]         = ax 
    dstates[7]         = dδf 
    
    return dstates
end

# cost weights
w_y  = 1e-5
w_sr = 2
w_ax = 0.2
w_ux = 0.2
w_sa = 1
function bicycle_cost(states, controls)
    y_cost = w_y * (states[2] - 5)^2
    dδf_cost = w_sr * (controls[2])^2
    ax_cost = w_ax * (controls[1])^2
    ux_cost = w_ux * (states[6] - 13.0)^2
    δf_cost = w_sa * (states[7])^2
    return y_cost + dδf_cost + ax_cost + ux_cost + δf_cost
end

using Plots

XL = [-40, -20, -3, -pi/5, -pi/2, 5.0, -pi/12]
XU = [300, 20, 3, pi/5, pi/2, 15.0, pi/12]
CL = [-2.6, -0.1]
CU = [2.6, 0.1]
X0 = [-10.0, 0, 0, 0, 0, 10.0, 0]
XF = [NaN, NaN, NaN, NaN, NaN, NaN, NaN]
ocp = defineOCP(numStates=7, numControls=2, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU);
defineStates!(ocp, [:x, :y, :v, :r, :ψ, :ux, :δf])
defineControls!(ocp, [:ax, :dδf])
OCPForm = ConfigurePredefined(ocp; (:Np => 81), (:tfDV => false), (:tf => 8), (:IntegrationScheme => :bkwEuler), (:dx => bicycleModel_expr), (:expr => bicycle_cost))
user_options = ()

OCPdef!(ocp, OCPForm)
xpos = ocp.p.x[:, 1];
y = ocp.p.x[:, 2]; dδf = ocp.p.u[:, 2]; ax = ocp.p.u[:, 1]; ux = ocp.p.x[:, 6]; δf = ocp.p.x[:, 7];
obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts], 36 <= ((xpos[i] - 30).^2 + (y[i] - 2).^2));
obj = @expression(ocp.f.mdl, sum((0.05 * (y[j] - sin(xpos[j]))^2 + 2 * dδf[j]^2 + 0.2 * ax[j]^2 + 0.2 * (ux[j] - 13)^2 + 1 * δf[j]^2) * ocp.f.TInt[j-1] for j in 2:ocp.f.Np))
@objective(ocp.f.mdl, Min, obj)
@time OptSolve!(ocp)
# plot(ocp.r.X[:, 1], ocp.r.X[:, 2], aspect_ratio = 1)
# # plot(ocp.r.X[:, 1], ocp.r.X[:, 3])

Dis = (ocp.r.X[:,1].-30).^2 + (ocp.r.X[:,2].-2).^2
@testset "bicycleTests" begin
    @test size(ocp.r.X)==(81, 7)
    @test count(x->x<35.9, Dis)==0
end