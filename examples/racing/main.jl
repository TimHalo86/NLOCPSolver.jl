include("../../src/utils.jl")
include("VehicleModel.jl")
include("racing_utils.jl")
using Plots
using MAT
using Interpolations
import Ipopt
import HSL_jll
using LinearAlgebra
using Statistics




function DefineOCP(cur_states, block_num, nCol, T)
    println(pwd())
    mat = matread("./processed_ThunderHill_West.mat");
    global center_line = transpose(mat["fined_center_line"])
    global lane_yaw = transpose(mat["fined_yawang"])
    left_lane_bound = transpose(mat["fined_left_bound"])
    right_lane_bound = transpose(mat["fined_right_bound"])
    block_info = transpose(mat["block_info"])


    ########################   Define Model Parameters ######################## 
    la      = 1.36
    lb      = 1.52
    hcg     = 0.525
    muf     = 0.98 # front friction
    mur     = 1.03 # rear friction
    M       = 2048 # Vehicle Mass
    g       = 9.81 # Gravity
    FZF0    = M * lb / (la + lb) * g
    FZR0    = M * la / (la + lb) * g
    KZF     = (hcg / (la + lb)) * M
    ratio   = 0.75
    ax_max  = (mur*FZR0 )/(M - mur*KZF) * 0.95
    ax_min  = max((-mur*FZR0 )/(M*(1-ratio) + mur*KZF), (-muf*FZF0 )/(M* ratio - muf*KZF)) * 0.95
    Caf     = 150000 # Front Cornering Stiffness
    Car     = 410000 # Rear Cornering Stiffness
    Izz     = 3675 # Vehicle Moment of inertia
    p = 4.0
    rhoKS = 5.0
    ######################## Initial States Settings ######################## 
    # states: x_s, y_s, v_s, r_s, ψ_s, ux_s, sa_s, ax_s, sr_s, jx_s, Δt_s

    ######################## Prepare cost_to_go Parameters and Block for KS based on current states ######################## 

    ######################## IMPORTANT!!!! Setting for IPOPT ######################## 
    user_options = (
        "mu_strategy" => "adaptive",
        # "linear_solver" => "ma27",
        "max_iter" => 150,
        "tol" => 4e-2,
        "dual_inf_tol" => 2.,
        "constr_viol_tol" => 5e-1,
        "compl_inf_tol" => 5e-1,
        "acceptable_tol" => 1.5e-1,
        "acceptable_constr_viol_tol" => 0.02,
        "acceptable_dual_inf_tol" => 1e10,
        "acceptable_compl_inf_tol" => 0.02,
        "warm_start_init_point" => "yes",
        "fixed_variable_treatment" => "relax_bounds",
        "max_cpu_time" => 0.1,
        "print_level" => 0,
    )

    XL = [NaN, NaN, -5, -2*pi, NaN, 0.1, -pi/3, ax_min]
    XU = [NaN, NaN, 5, 2*pi, NaN, 40, pi/3, ax_max]
    CL=[-0.15, -6]
    CU=[ 0.15, 6]

    init_pos = cur_states[1:2]'
    mpc_center_line = FindCenterLine(init_pos, center_line, lane_yaw)
    poly33_params = GetLengthParams(mpc_center_line)
    bg_idx, distance =  FindClosestPoint(init_pos, block_info[:,1:2])
    bg_idx = bg_idx - 2
    blocks = block_info[bg_idx:(bg_idx + block_num -1),:]
    block_yaw_info = blocks[:,3]
    block_center_x_info = blocks[:,1]
    block_center_y_info = blocks[:,2]
    block_length_info = blocks[:,4]
    block_width_info = blocks[:,5]
    block_s = zeros(block_num, 5)
    for num_idx = 1:1:block_num
        block_s[num_idx, 1] = block_center_x_info[num_idx]
        block_s[num_idx, 2] = block_center_y_info[num_idx]
        block_s[num_idx, 3] = block_yaw_info[num_idx]
        block_s[num_idx, 4] = block_width_info[num_idx]
        block_s[num_idx, 5] = block_length_info[num_idx]
    end

    ocp = defineOCP(numStates=8,numControls=2,X0=cur_states,XF=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], XL=XL, XU=XU, CL=CL,CU=CU);
    defineStates!(ocp, [:x,:y,:v,:r,:psi,:ux,:sa,:ax])
    defineControls!(ocp, [:sr, :jx])
    OCPForm = ConfigurePredefined(ocp; (:Np=>nCol), (:tfDV => false), (:tf => T), (:IntegrationScheme=>:bkwEuler), (:dx => LC500Model))
    # scheme = [repeat([:RK2], 2); repeat([:bkwEuler], nCol - 2)]
    # OCPForm.IntegrationScheme = scheme

    # nCol1 = 8
    # OCPForm.TInt = [0.1*ones(nCol1); ((T - 0.1*nCol1)/(nCol-nCol1-1))*ones(nCol-nCol1-1)]
    # OCPForm.tw = OCPForm.TInt./T

    OCPdef!(ocp, OCPForm)
    set_attributes(ocp.f.mdl, user_options...)

    x = ocp.p.x[:, 1]; 
    y = ocp.p.x[:, 2];
    v = ocp.p.x[:, 3]; 
    r = ocp.p.x[:, 4];
    ψ = ocp.p.x[:, 5]; 
    ux = ocp.p.x[:, 6];
    sa = ocp.p.x[:, 7]; 
    ax = ocp.p.x[:, 8];
    sr = ocp.p.u[:, 1]; 
    jx = ocp.p.u[:, 2];

    ocp.f.mdl[:block_center_x] = @variable(ocp.f.mdl, block_center_x[i = 1:block_num] in Parameter(block_center_x_info[i]))
    ocp.f.mdl[:block_center_y] = @variable(ocp.f.mdl, block_center_y[i = 1:block_num] in Parameter(block_center_y_info[i]))
    ocp.f.mdl[:block_yaw] = @variable(ocp.f.mdl, block_yaw[i = 1:block_num] in Parameter(block_yaw_info[i]))
    ocp.f.mdl[:block_width] = @variable(ocp.f.mdl, block_width[i = 1:block_num] in Parameter(block_width_info[i]))
    ocp.f.mdl[:block_length] = @variable(ocp.f.mdl, block_length[i = 1:block_num] in Parameter(block_length_info[i]))
    ocp.f.mdl[:cost_para] = @variable(ocp.f.mdl, cost_para[i = 1:10] in Parameter(poly33_params[i]))

    k_cost = @expression( ocp.f.mdl, sum((r[j]/ux[j])^2  for j=1:1:ocp.s.states.pts))
    v_cost = @expression( ocp.f.mdl, sum((v[j])^2  for j=1:1:ocp.s.states.pts))

    sr_cost = @expression( ocp.f.mdl, sum((sr[j])^2  for j=1:1:ocp.s.states.pts))
    jx_cost = @expression( ocp.f.mdl, sum((jx[j])^2  for j=1:1:ocp.s.states.pts))
    ax_cost = @expression( ocp.f.mdl, sum((ax[j])^2  for j=1:1:ocp.s.states.pts))
    sa_cost = @expression( ocp.f.mdl, sum((sa[j])^2  for j=1:1:ocp.s.states.pts))
    obj_goal = @expression(ocp.f.mdl, cost_para[1] + cost_para[2] * x[end] + cost_para[3] * y[end] + cost_para[4] * x[end]^2 + cost_para[5] * x[end] * y[end] + cost_para[6] * y[end]^2 + cost_para[7] * x[end]^3 + cost_para[8] * x[end]^2 * y[end] + cost_para[9] * x[end] * y[end]^2 + cost_para[10] * y[end]^3)

    obs_cost = @expression( ocp.f.mdl, sum(100*1/10*log(1 + exp(-10*(0.3 +(1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]  )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i] - 0.35) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )))  )) for j = 1:ocp.s.states.pts) )

    # MPC Original Constraint
    # DrivableTubeHard = @constraint(ocp.f.mdl, [j = 2:ocp.s.states.pts], 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )>= -0.2)

    # @constraint(ocp.f.mdl, [j = 2:ocp.s.states.pts], ax[j]  <= -0.1292 * (ux[j] - 60))

    # Discrete CBF Constraint
    # DrivableTubeHard = @constraint(ocp.f.mdl, [j = 1:ocp.s.states.pts-1], 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j+1]-block_center_x[i]) + sin(block_yaw[i]) * (y[j+1]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j+1]-block_center_x[i]) + cos(block_yaw[i]) * (y[j+1]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )
    #  - 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) ) 
    #  >= -0.5*(0.2 + 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )))

    # @constraint(ocp.f.mdl, [j = 1:ocp.s.states.pts-1], (-0.1292 * (ux[j+1] - 60) - ax[j+1])-(-0.1292 * (ux[j] - 60) - ax[j]) >= -0.5*(-0.1292 * (ux[j] - 60) - ax[j]))

    # Discrete CBF Constraint with a Compensator
    # alphaf = (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])
    # alphar = (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))
    # FYFMAX_d = (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))
    # FYRMAX_d = (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))
    # FYF_d = (-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))
    # FYR_d = (-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5))
    # fd = sqrt(((FYF+FYR)/M)^2+((la*FYF-lb*FYR)/Izz)^2)
    # fd = sqrt((((-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))+(-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5)))/M)^2+((la*(-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))-lb*(-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5)))/Izz)^2)
    DrivableTubeHard = @constraint(ocp.f.mdl, [j = 2:ocp.s.states.pts-1], 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j+1]-block_center_x[i]) + sin(block_yaw[i]) * (y[j+1]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j+1]-block_center_x[i]) + cos(block_yaw[i]) * (y[j+1]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )
    - 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) ) 
    >= -0.8*(0.2 + 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) ))
    + sqrt((((-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))+(-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5)))/M)^2+((la*(-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))-lb*(-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5)))/Izz)^2)
     * abs(1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j+1]-block_center_x[i]) + sin(block_yaw[i]) * (y[j+1]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j+1]-block_center_x[i]) + cos(block_yaw[i]) * (y[j+1]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )
    - 1/rhoKS *log(sum(exp( rhoKS*( - ( ((   (cos(block_yaw[i]) * (x[j]-block_center_x[i]) + sin(block_yaw[i]) * (y[j]-block_center_y[i]))/block_length[i]   )^p + ( (-sin(block_yaw[i]) * (x[j]-block_center_x[i]) + cos(block_yaw[i]) * (y[j]-block_center_y[i]))/(block_width[i]) )^p + 0.01)^(1/p) ) + 1 )  ) for i = 1:block_num) )))

    @constraint(ocp.f.mdl, [j = 2:ocp.s.states.pts-1], (-0.1292 * (ux[j+1] - 60) - ax[j+1])-(-0.1292 * (ux[j] - 60) - ax[j]) >= -0.8*(-0.1292 * (ux[j] - 60) - ax[j]) + abs((-0.1292 * (ux[j+1] - 60) - ax[j+1])-(-0.1292 * (ux[j] - 60) - ax[j])) * sqrt((((-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))+(-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5)))/M)^2+((la*(-2 * (sqrt((0.3 * (FZF0-KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Caf / (sqrt((0.3 * (FZF0-KZF*ax[j]))^2))) * (atan(( v[j] + la * r[j]) / (ux[j] + 0.1)) - sa[j])))-0.5))-lb*(-2 * (sqrt((0.3 * (FZR0+KZF*ax[j]))^2)) *(1/(1+exp(-(2 * Car / (sqrt((0.3 * (FZR0+KZF*ax[j]))^2))) * (atan((v[j] - lb * r[j]) / (ux[j] + 0.1)))))-0.5)))/Izz)^2))

    @objective(ocp.f.mdl, Min,  30*obs_cost + 12*obj_goal + 10*sr_cost + 10*sa_cost + 0.02*ax_cost + 0.01*jx_cost + 5*k_cost + 5*v_cost)
    # @time OptSolve!(ocp)
    return ocp
end


##################################################################
using Plots
plot()

mat = matread("processed_ThunderHill_West.mat");
center_line = transpose(mat["fined_center_line"])
lane_yaw = transpose(mat["fined_yawang"])
left_lane_bound = transpose(mat["fined_left_bound"])
right_lane_bound = transpose(mat["fined_right_bound"])
block_info = transpose(mat["block_info"])
block_num = 25 
nCol = 25

cur_states = [-1.82768, 50.1731, 0.0, 0.0, 1.6806, 5.0, 0.0, 0.0]
T = 4
ocp = DefineOCP(cur_states, block_num, nCol, T)


states_his = vcat(cur_states, 0, 0)
sol_t_list = []
t_sim = 0.0

δt_sim = 1e-2
states_real = cur_states[1:7]
states_real_his = states_real

OptSolve!(ocp)





for i = 1:1:2000
    println("time = ", t_sim)
    @time begin
        global center_line, lane_yaw, states_his, sol_t_list, t_sim, δt_sim, states_real, states_real_his, block_num,T
        x = ocp.p.x[:, 1]; 
        y = ocp.p.x[:, 2];
        v = ocp.p.x[:, 3]; 
        r = ocp.p.x[:, 4];
        ψ = ocp.p.x[:, 5]; 
        ux = ocp.p.x[:, 6];
        sa = ocp.p.x[:, 7]; 
        ax = ocp.p.x[:, 8];
        sr = ocp.p.u[:, 1]; 
        jx = ocp.p.u[:, 2];
        
        time_list = ocp.r.Tst
        time_list[end] = T
        time_list = time_list.+t_sim


        InterpolateAx = interpolate((time_list ,), value.(ax), Gridded(Constant{Next}()))
        InterpolateSr = interpolate((time_list ,), value.(sr), Gridded(Constant{Next}()))
        InterpolateJx = interpolate((time_list ,), value.(jx), Gridded(Constant{Next}()))

        global initial_guess = value.([x y v r ψ ux sa ax sr jx])
        cur_states = initial_guess[2,:]
        states_his = [states_his cur_states]

        for sim_idx in 1:Int32(floor(0.1/δt_sim))
            dstates = RacingVehicleModelPropagateTRI(states_real, [InterpolateSr(t_sim), InterpolateAx(t_sim)])
            states_real = states_real + dstates*δt_sim
            t_sim = t_sim + δt_sim
        end
        states_real_his = [states_real_his states_real]


        initial_guess[1:end-1,:] = initial_guess[2:end,:]
        initial_guess[1, 1:7] = states_real
        initial_guess[1, 8] = InterpolateAx(t_sim)
        initial_guess[1, 9] = InterpolateSr(t_sim)
        initial_guess[1, 10] = InterpolateJx(t_sim)
        
        mpc_center_line = FindCenterLine([states_real[1]; states_real[2]]', center_line, lane_yaw)
        poly33_params = GetLengthParams(mpc_center_line)
        bg_idx, distance = FindClosestPoint([states_real[1]; states_real[2]]', block_info[:,1:2])
        blocks = block_info[bg_idx:(bg_idx + block_num -1),:]
        block_yaw_info = blocks[:,3]
        block_center_x_info = blocks[:,1]
        block_center_y_info = blocks[:,2]
        block_length_info = blocks[:,4]
        block_width_info = blocks[:,5]

        block_s = zeros(block_num, 5)
        for num_idx = 1:1:block_num
            block_s[num_idx, 1] = block_center_x_info[num_idx]
            block_s[num_idx, 2] = block_center_y_info[num_idx]
            block_s[num_idx, 3] = block_yaw_info[num_idx]
            block_s[num_idx, 4] = block_width_info[num_idx]
            block_s[num_idx, 5] = block_length_info[num_idx]
        end

        set_parameter_value.(ocp.f.mdl[:block_center_x], block_center_x_info)
        set_parameter_value.(ocp.f.mdl[:block_center_y], block_center_y_info)
        set_parameter_value.(ocp.f.mdl[:block_yaw], block_yaw_info)
        set_parameter_value.(ocp.f.mdl[:block_width], block_width_info)
        set_parameter_value.(ocp.f.mdl[:block_length], block_length_info)
        set_parameter_value.(ocp.f.mdl[:cost_para], poly33_params)

        X0 = [states_real[1:7];initial_guess[1, 8]]
        U0 = [initial_guess[1,9];initial_guess[1,10]]
        UpdateX0!(ocp, X0)

        OptSolve!(ocp)
        sol_t_list = [sol_t_list; solve_time(ocp.f.mdl)]
    end
    println("ax = ", value(ax[1]))

    h1 = plot(value.(x), value.(y),aspect_ratio=:equal, lc=:green, marker = :circle, layout = (1, 2),subplot=1)
    plot!(h1, left_lane_bound[bg_idx*10:(bg_idx*10+280), 1], left_lane_bound[bg_idx*10:(bg_idx*10+280), 2],aspect_ratio=:equal, lc=:black, layout = (1, 2),subplot=1)
    plot!(h1, right_lane_bound[bg_idx*10:(bg_idx*10+280), 1], right_lane_bound[bg_idx*10:(bg_idx*10+280), 2],aspect_ratio=:equal, lc=:black, layout = (1, 2),subplot=1)
    plot!(h1, mpc_center_line[:,1], mpc_center_line[:,2], lc=:red, title =  "ux = $(round(2.23*initial_guess[1,6]; digits = 2)) mph", layout = (1, 2),subplot=1, size = (1000, 1000))
    
    plot!(h1, states_real_his[1,:], states_real_his[2,:], aspect_ratio=:equal, lc=:green, layout = (1, 2),subplot=2)
    plot!(h1, left_lane_bound[:, 1], left_lane_bound[:, 2],aspect_ratio=:equal, lc=:black, layout = (1, 2),subplot=2)
    plot!(h1, right_lane_bound[:, 1], right_lane_bound[:, 2],aspect_ratio=:equal, lc=:black, layout = (1, 2),subplot=2, size = (1000, 1000))
    
    display(h1)

end