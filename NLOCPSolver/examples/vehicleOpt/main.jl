# include("../../src/utils.jl")
# include("../../src/NLOCPSolver.jl")
using NLOCPSolver
include("bicycleModel.jl")
include("parameters.jl")
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
obj = @expression(ocp.f.mdl, sum((0.05 * (y[j])^2 + 2 * dδf[j]^2 + 0.2 * ax[j]^2 + 0.2 * (ux[j] - 13)^2 + 1 * δf[j]^2) * ocp.f.TInt[j-1] for j in 2:ocp.f.Np))
@objective(ocp.f.mdl, Min, obj)
@time OptSolve!(ocp)
plot(ocp.r.X[:, 1], ocp.r.X[:, 2], aspect_ratio = 1)
# plot(ocp.r.X[:, 6])
