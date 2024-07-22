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
OCPForm = ConfigurePredefined(ocp; (:Np => 81), (:tfDV => false), (:tf => 8), (:IntegrationScheme => :RK1), (:dx => bicycleModel_expr), (:expr => bicycle_cost))
# user_options = ()


OCPdef!(ocp, OCPForm)
block_center_x_info = 30;
block_center_y_info = 2;
block_radius_info = 6;
ocp.f.mdl[:block_center_x] = @variable(ocp.f.mdl, block_center_x in Parameter(block_center_x_info));
ocp.f.mdl[:block_center_y] = @variable(ocp.f.mdl, block_center_y in Parameter(block_center_y_info));
ocp.f.mdl[:block_radius] = @variable(ocp.f.mdl, block_radius in Parameter(block_radius_info));
xpos = ocp.p.x[:, 1];
y = ocp.p.x[:, 2]; dδf = ocp.p.u[:, 2]; ax = ocp.p.u[:, 1]; ux = ocp.p.x[:, 6]; δf = ocp.p.x[:, 7];
# obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts], 36 <= ((xpos[i] - 30).^2 + (y[i] - 2).^2));
set_parameter_value.(ocp.f.mdl[:block_center_x], block_center_x_info);
set_parameter_value.(ocp.f.mdl[:block_center_y], block_center_y_info);
set_parameter_value.(ocp.f.mdl[:block_radius], block_radius_info);
obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts],  block_radius^2 <= ((xpos[i] - block_center_x).^2 + (y[i] - block_center_y).^2));
# obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts-1], (xpos[i+1] - block_center_x).^2 + (y[i+1] - block_center_y).^2 - (xpos[i] - block_center_x).^2 - (y[i] - block_center_y).^2 >= -1 * ((xpos[i] - block_center_x).^2 + (y[i] - block_center_y).^2 - block_radius^2))
obj = @expression(ocp.f.mdl, sum((0.05 * (y[j])^2 + 2 * dδf[j]^2 + 0.2 * ax[j]^2 + 0.2 * (ux[j] - 13)^2 + 1 * δf[j]^2) * ocp.f.TInt[j-1] for j in 2:ocp.f.Np))
@objective(ocp.f.mdl, Min, obj)
@time OptSolve!(ocp)
function circleShape(h,k,r)
    theta = LinRange(0,2*pi,500)
    h .+ r*sin.(theta), k .+ r*cos.(theta)
end
plot()
plot!(circleShape(30,2,6), seriestype = [:shape,], lw=0.5, c=:blue, linecolor=:black, legend=false, fillalpha=0.5, aspect_ratio=1)
plot!(ocp.r.X[:, 1], ocp.r.X[:, 2], aspect_ratio = 1)
# plot(ocp.r.X[:, 6])
