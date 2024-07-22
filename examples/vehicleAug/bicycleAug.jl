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

    Fyf = MagicFormula(αf, Fzf, 1) # Front lateral force
    Fyr = MagicFormula(αr, Fzr, 1) # Rear lateral force

    dstates = Vector{Any}(undef,14)
    dstates[1]         = ux*cos(ψ) - v*sin(ψ) 
    dstates[2]         = ux*sin(ψ) + v*cos(ψ)  
    dstates[3]         = (Fyf+Fyr)/m - ux*r  
    dstates[4]         = (Fyf*la-Fyr*lb)/Izz 
    dstates[5]         = r
    dstates[6]         = ax 
    dstates[7]         = dδf 
    
    xa = states[8]
    ya = states[9]
    va = states[10]
    ra = states[11]
    ψa = states[12]
    uxa = states[13]
    δfa = states[14]
    axa = controls[3]
    dδfa = controls[4]
    Fzfa = m*g*la/(la+lb) - m*h/(la+lb)*axa # Front axle load
    Fzra = m*g*lb/(la+lb) + m*h/(la+lb)*axa # Rear axle load
    αfa = δfa - (va+la*ra)/uxa # Front slip angle
    αra = -(va-lb*ra)/uxa # Rear slip angle
    Fyfa = MagicFormula(αfa, Fzfa, 0.1) # Front lateral force
    Fyra = MagicFormula(αra, Fzra, 0.1) # Rear lateral force
    dstates[8]         = uxa*cos(ψa) - va*sin(ψa) 
    dstates[9]         = uxa*sin(ψa) + va*cos(ψa)  
    dstates[10]         = (Fyfa+Fyra)/m - uxa*ra  
    dstates[11]         = (Fyfa*la-Fyra*lb)/Izz
    dstates[12]         = ra 
    dstates[13]         = axa 
    dstates[14]         = dδfa 

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