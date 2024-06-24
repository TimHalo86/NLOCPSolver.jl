function ThreeDOFBicycle_expr(states, controls, parameter)
  # lateral tire load
  x = states[1]
  y = states[2]
  v = states[3]
  r = states[4]
  psi = states[5]
  sa = states[6]
  ux = states[7]
  ax = states[8]
  sr = controls[1]
  jx = controls[2]
  mu1 = parameter[1]
  mu2 = parameter[2]
  # mu = paramter[1]

  FYF=  ((PD2*(FzF0 - (ax - v*r)*KZX)^2 + PD1*(FzF0 - (ax - v*r)*KZX))*sin(PC1*atan((((PK1*sin(2*atan(PK2*(FzF0 - (ax - v*r)*KZX))))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)) + ((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX))^2 + EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)) - ((PE2*(FzF0 - (ax - v*r)*KZX) + PE1)*(1 - PE3)*(((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1))/((((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)^2 + EP_SMALL^2)^(0.5)))*((((PK1*sin(2*atan(PK2*(FzF0 - (ax - v*r)*KZX))))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)) + ((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX))^2 + EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)) - atan((((PK1*sin(2*atan(PK2*(FzF0 - (ax - v*r)*KZX))))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)) + ((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX))^2 + EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)))))) + (PV2*(FzF0 - (ax - v*r)*KZX)^2 + PV1*(FzF0 - (ax - v*r)*KZX)));
  FYR=  ((PD2*(FzR0 + (ax - v*r)*KZX)^2 + PD1*(FzR0 + (ax - v*r)*KZX))*sin(PC1*atan((((PK1*sin(2*atan(PK2*(FzR0 + (ax - v*r)*KZX))))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)) + ((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX))^2+EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)) - ((PE2*(FzR0 + (ax - v*r)*KZX) + PE1)*(1 - PE3*(((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1))/((((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)^2 + EP_SMALL^2)^(0.5))))*((((PK1*sin(2*atan(PK2*(FzR0 + (ax - v*r)*KZX))))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)) + ((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX))^2+EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)) - atan((((PK1*sin(2*atan(PK2*(FzR0 + (ax - v*r)*KZX))))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)) + ((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX))^2+EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)))))) + (PV2*(FzR0 + (ax - v*r)*KZX)^2 + PV1*(FzR0 + (ax - v*r)*KZX)));
  # FYF= mu1 * ((PD2*(FzF0 - (ax - v*r)*KZX)^2 + PD1*(FzF0 - (ax - v*r)*KZX))*sin(PC1*atan((((PK1*sin(2*atan(PK2*(FzF0 - (ax - v*r)*KZX))))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)) + ((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX))^2 + EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)) - ((PE2*(FzF0 - (ax - v*r)*KZX) + PE1)*(1 - PE3)*(((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1))/((((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)^2 + EP_SMALL^2)^(0.5)))*((((PK1*sin(2*atan(PK2*(FzF0 - (ax - v*r)*KZX))))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)) + ((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX))^2 + EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)) - atan((((PK1*sin(2*atan(PK2*(FzF0 - (ax - v*r)*KZX))))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)) + ((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX)))/(((PD2*PC1*(FzF0 - (ax - v*r)*KZX)^2 + PD1*PC1*(FzF0 - (ax - v*r)*KZX))^2 + EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v + la*r)/(ux+EP_SMALL)) - sa) + PH2*(FzF0 - (ax - v*r)*KZX) + PH1)))))) + (PV2*(FzF0 - (ax - v*r)*KZX)^2 + PV1*(FzF0 - (ax - v*r)*KZX)));
  # FYR= mu2 * ((PD2*(FzR0 + (ax - v*r)*KZX)^2 + PD1*(FzR0 + (ax - v*r)*KZX))*sin(PC1*atan((((PK1*sin(2*atan(PK2*(FzR0 + (ax - v*r)*KZX))))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)) + ((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX))^2+EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)) - ((PE2*(FzR0 + (ax - v*r)*KZX) + PE1)*(1 - PE3*(((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1))/((((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)^2 + EP_SMALL^2)^(0.5))))*((((PK1*sin(2*atan(PK2*(FzR0 + (ax - v*r)*KZX))))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)) + ((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX))^2+EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)) - atan((((PK1*sin(2*atan(PK2*(FzR0 + (ax - v*r)*KZX))))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)) + ((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX)))/(((PD2*PC1*(FzR0 + (ax - v*r)*KZX)^2 + PD1*PC1*(FzR0 + (ax - v*r)*KZX))^2+EP_SMALL^2)^(0.5))*0.001)+EP_SMALL))*((atan((v - lb*r)/(ux+EP_SMALL))) + PH2*(FzR0 + (ax - v*r)*KZX) + PH1)))))) + (PV2*(FzR0 + (ax - v*r)*KZX)^2 + PV1*(FzR0 + (ax - v*r)*KZX)));

  dx = Vector{Any}(undef,8)
  dx[1] = (ux*cos(psi) - (v + la*r)*sin(psi))   # X position
  dx[2] = (ux*sin(psi) + (v + la*r)*cos(psi))   # Y position
  dx[3] = ((FYF + FYR)/m - r*ux)                       # Lateral Speed
  dx[4] = ((la*FYF-lb*FYR)/Izz)                            # Yaw Rate
  dx[5] = (r)                                                # Yaw Angle
  dx[6] = (sr)                                               # Steering Angle
  dx[7] = (ax)                                               # Longitudinal Speed
  dx[8] = (jx)                                               # Longitudinal Acceleration
  return dx
end


function ThreeDOFBicycle_cons(states, controls, parameter)
  ux = states[7]
  ax = states[8]
  cons = [ax + 0.1296 * (ux - 60); ax + 0.1296 * (ux - 50)]
  return cons
end

function ThreeDOFBicycle_cost(states, controls, parameter)
  sr = controls[2]
  return sr^2
end