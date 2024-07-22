# include("../../src/utils.jl")
# include("../../src/NLOCPSolver.jl")
using NLOCPSolver
include("bicycleAug.jl")
using Plots

XL = [-40, -20, -3, -pi/5, -pi/2, 5.0, -pi/12, -40, -20, -3, -pi/5, -pi/2, 5.0, -pi/12]
XU = [300, 20, 3, pi/5, pi/2, 15.0, pi/12, 300, 20, 3, pi/5, pi/2, 15.0, pi/12]
CL = [-2.6, -0.1, -2.6, -0.1]
CU = [2.6, 0.1, 2.6, 0.1]
X0 = [-10.0, 0, 0, 0, 0, 10.0, 0, -10.0, 0, 0, 0, 0, 10.0, 0]
XF = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
ocp = defineOCP(numStates=14, numControls=4, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU);
defineStates!(ocp, [:x, :y, :v, :r, :ψ, :ux, :δf, :xa, :ya, :va, :ra, :ψa, :uxa, :δfa])
defineControls!(ocp, [:ax, :dδf, :axa, :dδfa])
OCPForm = ConfigurePredefined(ocp; (:Np => 81), (:tfDV => false), (:tf => 8), (:IntegrationScheme => :bkwEuler), (:dx => bicycleModel_expr), (:expr => bicycle_cost))


OCPdef!(ocp, OCPForm)
xpos = ocp.p.x[:, 1];
y = ocp.p.x[:, 2]; dδf = ocp.p.u[:, 2]; ax = ocp.p.u[:, 1]; ux = ocp.p.x[:, 6]; δf = ocp.p.x[:, 7];
xposa = ocp.p.x[:, 8];
ya = ocp.p.x[:, 9]; dδfa = ocp.p.u[:, 4]; axa = ocp.p.u[:, 3]; uxa = ocp.p.x[:, 13]; δfa = ocp.p.x[:, 14];
# obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts], 36 <= ((xpos[i] - 30).^2 + (y[i] - 2).^2));
@constraint(ocp.f.mdl, [i=1:ocp.s.states.pts-1], (xpos[i+1] - 30).^2 + (y[i+1] - 2).^2 - (xpos[i] - 30).^2 - (y[i] - 2).^2 >= -1 * ((xpos[i] - 30).^2 + (y[i] - 2).^2 - 36))
@constraint(ocp.f.mdl, [i=1:ocp.s.states.pts-1], (xposa[i+1] - 30).^2 + (ya[i+1] - 2).^2 - (xposa[i] - 30).^2 - (ya[i] - 2).^2 >= -1 * ((xposa[i] - 30).^2 + (ya[i] - 2).^2 - 36))

obj = @expression(ocp.f.mdl, sum((0.05 * ((y[j])^2 + (ya[j])^2) + 2 * (dδf[j]^2 + dδfa[j]^2) + 0.2 * (ax[j]^2 + axa[j]^2) + 0.2 * ((ux[j] - 13)^2 + (uxa[j] - 13)^2) + 1 * (δf[j]^2 + δfa[j]^2)) * ocp.f.TInt[j-1] for j in 2:ocp.f.Np))
@objective(ocp.f.mdl, Min, obj)
@time OptSolve!(ocp)
function circleShape(h,k,r)
    theta = LinRange(0,2*pi,500)
    h .+ r*sin.(theta), k .+ r*cos.(theta)
end
# plot(1)
# plot!(ocp.r.X[:, 1], ocp.r.X[:, 2], label="mu=1", aspect_ratio = 1)
# plot!(ocp.r.X[:, 8], ocp.r.X[:, 9], label="mu=0.1", aspect_ratio = 1)
# plot!(circleShape(30,2,6), seriestype = [:shape,], lw=0.5, c=:blue, linecolor=:black, legend=false, fillalpha=0.5, aspect_ratio=1)

la = 1.56 
lb = 1.64
m = 2020
g = 9.81
Izz = 4095
h = 0.6
mu = 0.5
B = 5.68; C = 1.817;
U_sim = 1/2*(ocp.r.U[:, 1:2] + ocp.r.U[:, 3:4]);
ax_sim = U_sim[:, 1]; dδf_sim = U_sim[:, 2];
x0_sim = [-10.0, 0, 0, 0, 0, 10.0, 0];
x_sim = zeros(1,81); x_sim[1] = -10.0;
y_sim = zeros(1,81); y_sim[1] = 0;
v_sim = zeros(1,81); v_sim[1] = 0;
r_sim = zeros(1,81); r_sim[1] = 0;
ψ_sim = zeros(1,81); ψ_sim[1] = 0;
ux_sim = zeros(1,81); ux_sim[1] = 10.0;
δf_sim = zeros(1,81); δf_sim[1] = 0;
for i in 1:80
    Fzf = m*g*la/(la+lb) - m*h/(la+lb)*ax_sim[i]; # Front axle load
    Fzr = m*g*lb/(la+lb) + m*h/(la+lb)*ax_sim[i]; # Rear axle load
    αf = δf_sim[i] - (v_sim[i]+la*r_sim[i])/ux_sim[i]; # Front slip angle
    αr = -(v_sim[i]-lb*r_sim[i])/ux_sim[i]; # Rear slip angle
    Fyf = mu*Fzf*sin(C*atan(B/mu*αf));
    Fyr = mu*Fzr*sin(C*atan(B/mu*αr));
    x_sim[i+1] = x_sim[i] + (ux_sim[i]*cos(ψ_sim[i]) - v_sim[i]*sin(ψ_sim[i]))*0.1; 
    y_sim[i+1] = y_sim[i] + (ux_sim[i]*sin(ψ_sim[i]) + v_sim[i]*cos(ψ_sim[i]))*0.1;
    v_sim[i+1] = v_sim[i] + ((Fyf+Fyr)/m - ux_sim[i]*r_sim[i])*0.1;
    r_sim[i+1] = r_sim[i] + ((Fyf*la-Fyr*lb)/Izz)*0.1;
    ψ_sim[i+1] = ψ_sim[i] + r_sim[i]*0.1;
    ux_sim[i+1] = ux_sim[i] + ax_sim[i]*0.1;
    δf_sim[i+1] = δf_sim[i] + dδf_sim[i]*0.1;
end
plot(2)
plot!(ocp.r.X[:, 1], ocp.r.X[:, 2], label="mu=1", aspect_ratio = 1)
plot!(ocp.r.X[:, 8], ocp.r.X[:, 9], label="mu=0.1", aspect_ratio = 1)
plot!(vec(x_sim), vec(y_sim), label="mu=0.5", aspect_ratio = 1)
plot!(circleShape(30,2,6), label=false, seriestype = [:shape,], lw=0.5, c=:blue, linecolor=:black, legend=false, fillalpha=0.5, aspect_ratio=1)
