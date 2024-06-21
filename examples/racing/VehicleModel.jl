function RacingVehicleModelPropagateTRIESP(states, ctrl)
    x = states[1]; y = states[2]; v = states[3]; r = states[4]; psi = states[5]; ux = states[6]; sa = states[7];
    sr = ctrl[1]; ax = ctrl[2]; 

    alpha_fl = (atan(( v + la * r) / (ux + 0.1)) - sa)
    alpha_fr = (atan(( v + la * r) / (ux + 0.1)) - sa)

    alpha_rl = (atan(v - lb * r) / (ux + 0.1))
    alpha_rr = (atan(v - lb * r) / (ux + 0.1))

    FZFL = (FZF0 - KZF * (ax)) / 2 
    FZFR = (FZF0 - KZF * (ax)) / 2 

    FZRL = (FZR0 + KZF * (ax)) / 2 
    FZRR = (FZR0 + KZF * (ax)) / 2 

    if ax >= 0
        axf = 0
    else
        axf = ax * ratio
    end

    # axf = ((1 - 1/(1+exp(-10 * ax))) * ax * ratio)
    axr = (ax - axf)
    FXF = (axf * M)
    FXFL = FXF * FZFL / (FZFL + FZFR)
    FXFR = FXF * FZFR / (FZFL + FZFR)
    FXR = (axr * M)
    FXRL = FXR * FZRL / (FZRL + FZRR)
    FXRR = FXR * FZRR / (FZRL + FZRR)


    rdes = ux * sa / (la + lb);
    ΔFXFL = 0
    ΔFXFR = 0 
    α = Izz*0.2


    FXFL = FXFL + ΔFXFL
    FXFR = FXFR + ΔFXFR
    FXRL = FXRL - (ΔFXFR+ΔFXFR)/2
    FXRR = FXRR - (ΔFXFR+ΔFXFR)/2



    FYFLMAX = (sqrt(max((muf * FZFL)^2 - FXFL^2,  0.1)))
    FYFRMAX = (sqrt(max((muf * FZFR)^2 - FXFR^2,  0.1)))
    FYRLMAX = (sqrt(max((mur * FZRL)^2 - FXRL^2,  0.1)))
    FYRRMAX = (sqrt(max((mur * FZRR)^2 - FXRR^2,  0.1)))

    
    FYFL = (-2 * FYFLMAX *(1/(1+exp(-(1 * Caf / FYFLMAX) * alpha_fl))-0.5))
    FYFR = (-2 * FYFRMAX *(1/(1+exp(-(1 * Caf / FYFRMAX) * alpha_fr))-0.5))
    FYRL = (-2 * FYRLMAX *(1/(1+exp(-(1 * Car / FYRLMAX) * alpha_rl))-0.5))
    FYRR = (-2 * FYRRMAX *(1/(1+exp(-(1 * Car / FYRRMAX) * alpha_rr))-0.5))

    FX1         = FXFL*cos(sa) - FYFL*sin(sa); 
    FX2         = FXFR*cos(sa) - FYFR*sin(sa);
    FX3         = FXRL;
    FX4         = FXRR;
    FY1         = FXFL*sin(sa) + FYFL*cos(sa);
    FY2         = FXFR*sin(sa) + FYFR*cos(sa);
    FY3         = FYRL;
    FY4         = FYRR;

    dx = Array{Float64}(undef,7)
    dx[1] = (ux * cos(psi) - v * sin(psi))
    dx[2] = (ux * sin(psi) + v * cos(psi))
    dx[3] = (FY1 + FY2 + FY3 + FY4 ) / M -  r * ux
    dx[4] = ((FY1 + FY2) * la - (FY3 + FY4) * lb - (FX1 + FX3 - FX2 - FX4) * lt / 2) / Izz
    dx[5] = (r)
    dx[6] = (ax + v * r)
    dx[7] = (sr)
    return dx
end

function RacingVehicleModelPropagateTRI(states, ctrl)
    la      = 1.36
    lb      = 1.52
    hcg     = 0.525
    Izz     = 3675 # Vehicle Moment of inertia
    muf     = 0.98 # front friction
    mur     = 1.03 # rear friction
    M       = 2048 # Vehicle Mass
    g       = 9.81 # Gravity
    Caf     = 150000 # Front Cornering Stiffness
    Car     = 410000 # Rear Cornering Stiffness
    FZF0    = M * lb / (la + lb) * g
    FZR0    = M * la / (la + lb) * g
    KZF     = (hcg / (la + lb)) * M
    ratio   = 0.75

    
    
    x = states[1]; y = states[2]; v = states[3]; r = states[4]; psi = states[5]; ux = states[6]; sa = states[7];
    ax = ctrl[2]; sr = ctrl[1];
    alphaf = (atan(( v + la * r) / (ux + 0.1)) - sa)
    alphar = (atan((v - lb * r) / (ux + 0.1)))

    FZF = (FZF0 - KZF * (ax))
    FZR = (FZR0 + KZF * (ax))
    if ax >= 0
        axf = 0
    else
        axf = ax * ratio
    end
    axr = (ax - axf)
    FXF = (axf * M)
    FXR = (axr * M)
    FYFMAX = (sqrt(max((muf * FZF)^2 - FXF^2,  0.1)))
    FYRMAX = (sqrt(max((mur * FZR)^2 - FXR^2,  0.1)))
    
    FYF = (-2 * FYFMAX *(1/(1+exp(-(2 * Caf / FYFMAX) * alphaf))-0.5))
    FYR = (-2 * FYRMAX *(1/(1+exp(-(2 * Car / FYRMAX) * alphar))-0.5))
    dx = Array{Float64}(undef,7)
    dx[1] = (ux * cos(psi) - v * sin(psi))
    dx[2] = (ux * sin(psi) + v * cos(psi))
    dx[3] = ((FYF + FYR)/M - r*ux)
    dx[4] = (((la * FYF- lb * FYR)/Izz))
    dx[5] = (r)
    dx[6] = (ax + v * r)
    dx[7] = (sr)
    return dx
end


function LC500Model(states, ctrls, params)
    la      = 1.36
    lb      = 1.52
    hcg     = 0.525
    Izz     = 3675 # Vehicle Moment of inertia
    muf     = 0.98 # front friction
    mur     = 1.03 # rear friction
    M       = 2048 # Vehicle Mass
    g       = 9.81 # Gravity
    Caf     = 150000 # Front Cornering Stiffness
    Car     = 410000 # Rear Cornering Stiffness
    FZF0    = M * lb / (la + lb) * g
    FZR0    = M * la / (la + lb) * g
    KZF     = (hcg / (la + lb)) * M
    ratio   = 0.75

    v = states[3]; r = states[4]; psi = states[5]; ux = states[6]; sa = states[7]; ax = states[8];
    sr = ctrls[1]; jx = ctrls[2]

    alphaf = (atan(( v + la * r) / (ux + 0.1)) - sa)
    alphar = (atan(v - lb * r) / (ux + 0.1))

    
    axf = ((1 - 1/(1+exp(-10 * ax))) * ax * ratio) #(ax * ratio)
    axr = (ax - axf)

    FZF = (FZF0 - KZF * (ax ))
    FZR = (FZR0 + KZF * (ax ))
    FXF = (axf* M)
    FXR = (axr* M)

    FYFMAX = (sqrt( log(1+exp( 10.0*(1 - (FXF/(muf*FZF))^2 - 0.06) ))/10.0)*(muf*FZF)+ 0.1) 
    FYRMAX = (sqrt( log(1+exp( 10.0*(1 - (FXR/(mur*FZR))^2  - 0.06)  ))/10.0)*(mur*FZR)+ 0.1) 

    FYF = (-2 * FYFMAX *(1/(1+exp(-(2 * Caf / FYFMAX) * alphaf))-0.5))
    FYR = (-2 * FYRMAX *(1/(1+exp(-(2 * Car / FYRMAX) * alphar))-0.5))

    dx = Vector{Any}(undef,8)
    dx[1] = (ux * cos(psi) - v * sin(psi))
    dx[2] = (ux * sin(psi) + v * cos(psi))
    dx[3] = ((FYF *  cos(sa) + FYR)/M - r*ux)
    dx[4] = (((la * FYF * cos(sa) - lb * FYR)/Izz))
    dx[5] = (r)
    dx[6] = (ax+ r*v)
    dx[7] = (sr)
    dx[8] = (jx)
    return dx
end


# function ThreeDOFBicycle_cons(states, controls)
#     ux = states[6]; ax = states[8];
#     cons = [ax  <= -0.1292 * (ux - 60)]
#     return cons
# end


# function ThreeDOFBicycle_cost(states, ctrls)
#     v = states[3]; r = states[4]; psi = states[5]; ux = states[6]; sa = states[7]; ax = states[8];
#     sr = ctrls[1]; jx = ctrls[2]
#     return 10*sr^2 + 10*sa^2 + 0.02*ax^2 + 0.01*jx^2 + 5*(r/ux)^2 +5*v^2 
#   end