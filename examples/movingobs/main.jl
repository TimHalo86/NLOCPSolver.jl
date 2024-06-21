include("../../src/utils.jl")
include("parameter.jl")
include("ThreeDOF_Bicycle.jl")
using Plots


XL = [-2.6, 0, NaN, NaN, -6.283185307179586, -0.62, 10, -2.6]
XU = [1.801, 300, NaN, NaN, 6.283185307179586, 0.62, 20,  2.6]
CL=[-2.1, -6.2]
CU=[1.9, 6.2]
ocp = defineOCP(numStates=8,numControls=2,X0=[1.8, 0., 0.0, 0.0, pi/2, 0.0, 15., 0.],XF=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], XL=XL, XU=XU, CL=CL,CU=CU);
defineStates!(ocp, [:x,:y,:v,:r,:psi,:sa,:ux,:ax])
defineControls!(ocp, [:sr, :jx])
# trapezoidal
# bkwEuler
OCPForm = ConfigurePredefined(ocp; (:Np=>25), (:tfDV => true),  (:IntegrationScheme=>:RK2), (:dx => ThreeDOFBicycle_expr), (:expr=>ThreeDOFBicycle_cost), (:params => [0.3, 0.5]))
user_options = ()

# scheme = [repeat([:RK4], 5); repeat([:RK3], 5); repeat([:RK2], 5);repeat([:RK1], 5);repeat([:bkwEuler], 5);repeat([:trapezoidal], 5);]
# scheme = [repeat([:RK4], 3); repeat([:bkwEuler], 27)]
# scheme =  repeat([:bkwEuler], 30)
# OCPForm.IntegrationScheme = scheme


OCPdef!(ocp, OCPForm)
x = ocp.p.x[:, 1]; y = ocp.p.x[:, 2]; ux = ocp.p.x[:, 7]; sr = ocp.p.u[:, 1]; v = ocp.p.x[:, 3]
timeSeq = ocp.p.tV
obs_info = [1.8 50.7368 0.3 0.9 0.0 5.5]
ux_con = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts-1], (ux[i+1]-15)^2<=0.01)
lsm = 4.66 #to leave 2m gap between the tip of the vehicle and the behind of the bicycle
ssm = 3.3 #to move the vehicle to the center of the left lane.
obs_ind = 1
obs_con1 = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts-1], 1 <= ((x[(i+1)]-timeSeq[i+1]*obs_info[obs_ind, 5]-obs_info[obs_ind, 1])^2)/((obs_info[obs_ind, 3]+ssm)^2) + ((y[(i+1)]-timeSeq[i+1]*obs_info[obs_ind, 6]-obs_info[obs_ind, 2])^2)/((obs_info[obs_ind, 4]+lsm)^2));
obj = @expression(ocp.f.mdl,  sum((5 * (x[j] - 1.8)^2 + 10 * v[j]^2 +  10 * sr[j]^2 ) * ocp.f.TInt[j - 1] for  j in 2:ocp.f.Np) )
@objective(ocp.f.mdl, Min,  obj + ocp.f.tf + (y[end] - 200)^2 + (x[end] - 1.8)^2)
@time OptSolve!(ocp)
plot(ocp.r.X[:, 1], ocp.r.X[:, 2])