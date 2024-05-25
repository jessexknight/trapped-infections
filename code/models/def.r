library('deSolve')
# ------------------------------------------------------------------------------
# indices
# g = gen pop; k = key pop
# S = susceptible; A = acute infection; Y = chronic infection
# u = unpaired (all models); p = post-transm (EPA); gXX: pair states (pair-based)
ylabs = c(
  'gSu','gAu','gYu','kSu','kAu','kYu', # flat  1: 6
  'gAp','gYp','kAp','kYp',             # epa   7:10
  'gSS','gSA','gSY','gAA','gAY','gYY') # pair 11:16
nc = length(ylabs)
for (i in seq(nc)){ assign(ylabs[i],i,env=.GlobalEnv) } # HACK
# ------------------------------------------------------------------------------
# parameters
# unit: years
par.def = function(...){
  p = list(
    # populations
    prop_k = .02, # prop key pop entering
    size = 10000, # total size
    exit = 1/35, # pop turnover
    turn = 1/8.2, # group turnover
    # partnerships
    acts_sh = 3, # acts per short ptr
    rate_sh_kg = 60, # short ptrs rate: key with gen
    rate_sh_kk = 60, # short ptrs rate: key with key
    gap_lo_gg = .5, # mean gap b/w long ptrs
    dur_lo  = 10, # duration of long ptrs
    freq_lo = 70, # act freq per long ptr
    # infection
    beta_Y = .0007, # chronic transm risk per act
    beta_A = .0007*26, # acute transm risk per act
    mort = 1/11, # mort rate among chronic
    prog = 4, # acute-to-chronic rate
    init_Y = 1 # init num infections
  )
  # allow override
  pu = list(...)
  for (name in names(pu)){
    p[[name]] = pu[[name]]
  }
  # conditional
  p$prop_g = 1-p$prop_k # prop gen pop entering
  p$size_k = p$size * p$prop_k * p$exit / (p$turn + p$exit) # equilib key pop size
  p$size_g = p$size - p$size_k # equilib gen pop size
  p$rate_sh_gk = p$rate_sh_kg * p$size_k / p$size_g # short ptrs rate: gen with key
  p$rate_lo_gg = 1/(p$gap_lo_gg+p$dur_lo) # long ptrs rate: gen with gen
  p$cur_lo_gg = p$rate_lo_gg * p$dur_lo # mean num long ptrs among gen pop
  p$acts_lo = p$freq_lo * p$dur_lo # acts per long ptr
  p$prop_pair = 2/p$gap_lo_gg / (1/p$dur_lo + 2/p$gap_lo_gg + p$exit) # equilib prop paired
  # return
  return(p)
}
# ------------------------------------------------------------------------------
# flat
y0.flat = function(p){
  y0 = setNames(rep(0,nc),ylabs)
  y0[gSu] = p$size_g - p$init_Y
  y0[gYu] = p$init_Y
  y0[kSu] = p$size_k - p$init_Y
  y0[kYu] = p$init_Y
  return(y0)
}
dy.flat = function(t,y,p){
  dy = y*0
  # sums
  size_g = y[gSu]+y[gAu]+y[gYu]
  size_k = y[kSu]+y[kAu]+y[kYu]
  prev_g = (size_g-y[gSu]) / size_g
  prev_k = (size_k-y[kSu]) / size_k
  # foi
  B_A_sh = 1-(1-p$beta_A)^(p$acts_sh)
  B_Y_sh = 1-(1-p$beta_Y)^(p$acts_sh)
  B_A_lo = 1-(1-p$beta_A)^(p$acts_lo)
  B_Y_lo = 1-(1-p$beta_Y)^(p$acts_lo)
  foi_g = p$rate_lo_gg * (y[gAu]*B_A_lo + y[gYu]*B_Y_lo) / size_g +
          p$rate_sh_gk * (y[kAu]*B_A_sh + y[kYu]*B_Y_sh) / size_k
  foi_k = p$rate_sh_kg * (y[gAu]*B_A_sh + y[gYu]*B_Y_sh) / size_g +
          p$rate_sh_kk * (y[kAu]*B_A_sh + y[kYu]*B_Y_sh) / size_k
  # dy
  dy[gSu] = -y[gSu]*p$exit +y[kSu]*p$turn +p$size*p$prop_g*p$exit -y[gSu]*foi_g
  dy[kSu] = -y[kSu]*p$exit -y[kSu]*p$turn +p$size*p$prop_k*p$exit -y[kSu]*foi_k
  dy[gAu] = -y[gAu]*p$exit +y[kAu]*p$turn -y[gAu]*p$prog +y[gSu]*foi_g
  dy[kAu] = -y[kAu]*p$exit -y[kAu]*p$turn -y[kAu]*p$prog +y[kSu]*foi_k
  dy[gYu] = -y[gYu]*p$exit +y[kYu]*p$turn -y[gYu]*p$mort +y[gAu]*p$prog
  dy[kYu] = -y[kYu]*p$exit -y[kYu]*p$turn -y[kYu]*p$mort +y[kAu]*p$prog
  # return
  return(list(dy=dy,
    foi_g  = unname(foi_g),
    foi_k  = unname(foi_k),
    prev_g = unname(prev_g),
    prev_k = unname(prev_k)))
}
# ------------------------------------------------------------------------------
# epa
y0.epa = function(p){
  return(y0.flat(p))
}
dy.epa = function(t,y,p){
  dy = y*0
  # sums
  size_g = y[gSu]+y[gAu]+y[gAp]+y[gYu]+y[gYp]
  size_k = y[kSu]+y[kAu]+y[kAp]+y[kYu]+y[kYp]
  prev_g = (size_g-y[gSu]) / size_g
  prev_k = (size_k-y[kSu]) / size_k
  # foi: long (epa)
  Sgg_lo = y[gSu] * p$cur_lo_gg # cur eff ptrs among sus
  Agg_lo = y[gAu] * p$cur_lo_gg + y[gAp] * (p$cur_lo_gg-1) # cur eff ptrs among acute
  Ygg_lo = y[gYu] * p$cur_lo_gg + y[gYp] * (p$cur_lo_gg-1) # cur eff ptrs among chronic
  tot_lo = Sgg_lo + Agg_lo + Ygg_lo # total cur eff ptrs
  foi_Agg = Sgg_lo * Agg_lo / tot_lo * p$freq_lo * p$beta_A # abs trans from acute
  foi_Ygg = Sgg_lo * Ygg_lo / tot_lo * p$freq_lo * p$beta_Y # abs trans from chronic
  foi_Sgg = foi_Agg + foi_Ygg # abs trans to sus
  # foi: short (classic)
  B_A_sh = 1-(1-p$beta_A)^(p$acts_sh)
  B_Y_sh = 1-(1-p$beta_Y)^(p$acts_sh)
  foi_g = p$rate_sh_gk * ((y[kAu]+y[kAp])*B_A_sh + (y[kYu]+y[kYp])*B_Y_sh) / size_k +
          foi_Sgg / y[gSu]
  foi_k = p$rate_sh_kg * ((y[gAu]+y[gAp])*B_A_sh + (y[gYu]+y[gYp])*B_Y_sh) / size_g +
          p$rate_sh_kk * ((y[kAu]+y[kAp])*B_A_sh + (y[kYu]+y[kYp])*B_Y_sh) / size_k
  # dy
  dy[gSu] = -y[gSu]*p$exit +y[kSu]*p$turn +p$size*p$prop_g*p$exit         -y[gSu]*foi_g
  dy[kSu] = -y[kSu]*p$exit -y[kSu]*p$turn +p$size*p$prop_k*p$exit         -y[kSu]*foi_k
  dy[gAp] = -y[gAp]*p$exit +y[kAp]*p$turn -y[gAp]/p$dur_lo -y[gAp]*p$prog +y[gSu]*foi_g  +foi_Agg
  dy[kAp] = -y[kAp]*p$exit -y[kAp]*p$turn -y[kAp]/p$dur_lo -y[kAp]*p$prog +y[kSu]*foi_k
  dy[gAu] = -y[gAu]*p$exit +y[kAu]*p$turn +y[gAp]/p$dur_lo -y[gAu]*p$prog                -foi_Agg
  dy[kAu] = -y[kAu]*p$exit -y[kAu]*p$turn +y[kAp]/p$dur_lo -y[kAu]*p$prog
  dy[gYp] = -y[gYp]*p$exit +y[kYp]*p$turn -y[gYp]/p$dur_lo +y[gAp]*p$prog -y[gYp]*p$mort +foi_Ygg
  dy[kYp] = -y[kYp]*p$exit -y[kYp]*p$turn -y[kYp]/p$dur_lo +y[kAp]*p$prog -y[kYp]*p$mort
  dy[gYu] = -y[gYu]*p$exit +y[kYu]*p$turn +y[gYp]/p$dur_lo +y[gAu]*p$prog -y[gYu]*p$mort -foi_Ygg
  dy[kYu] = -y[kYu]*p$exit -y[kYu]*p$turn +y[kYp]/p$dur_lo +y[kAu]*p$prog -y[kYu]*p$mort
  # return
  return(list(dy=dy,
    foi_g  = unname(foi_g),
    foi_k  = unname(foi_k),
    prev_g = unname(prev_g),
    prev_k = unname(prev_k)))
}
# ------------------------------------------------------------------------------
# pair
y0.pair = function(p){
  y0 = setNames(rep(0,nc),ylabs)
  y0[gSu] = p$size_g * (1-p$prop_pair)
  y0[gSS] = p$size_g * p$prop_pair * 0.5 - p$init_Y
  y0[gSY] = p$init_Y
  y0[kSu] = p$size_k - p$init_Y
  y0[kYu] = p$init_Y
  return(y0)
}
dy.pair = function(t,y,p){
  dy = y*0
  # sums
  size_k  = y[kSu]+y[kAu]+y[kYu]
  size_gu = y[gSu]+y[gAu]+y[gYu]+1e-14
  size_gp = 2*(y[gSS]+y[gSA]+y[gSY]+y[gAA]+y[gAY]+y[gYY])+1e-14
  size_gS = y[gSu] + y[gSS]*2 + y[gSA] + y[gSY]
  size_g  = size_gu + size_gp
  prev_g  = (size_g-size_gS)/size_g
  prev_k  = (size_k-y[kSu])/size_k
  # form & end
  form_gu = 1/p$gap_lo_gg/size_gu
  form_SS = y[gSu]*y[gSu]*form_gu
  form_SA = y[gSu]*y[gAu]*form_gu
  form_SY = y[gSu]*y[gYu]*form_gu
  form_AA = y[gAu]*y[gAu]*form_gu
  form_AY = y[gAu]*y[gYu]*form_gu
  form_YY = y[gYu]*y[gYu]*form_gu
  form_Sx = form_SS*2 +form_SA   +form_SY
  form_Ax = form_SA   +form_AA*2 +form_AY
  form_Yx = form_SY   +form_AY   +form_YY*2
  end_Sx  = (y[gSS]*2 +y[gSA]   +y[gSY]  )/p$dur_lo
  end_Ax  = (y[gSA]   +y[gAA]*2 +y[gAY]  )/p$dur_lo
  end_Yx  = (y[gSY]   +y[gAY]   +y[gYY]*2)/p$dur_lo
  # foi
  B_A_sh = 1-(1-p$beta_A)^(p$acts_sh)
  B_Y_sh = 1-(1-p$beta_Y)^(p$acts_sh)
  y_gA = y[gAu] + y[gSA] + y[gAA]*2 + y[gAY]
  y_gY = y[gYu] + y[gSY] + y[gAY] + y[gYY]*2
  foi_Agg = p$freq_lo * p$beta_A * y[gSA]
  foi_Ygg = p$freq_lo * p$beta_Y * y[gSY]
  foi_gk = p$rate_sh_gk * (y[kAu]*B_A_sh + y[kYu]*B_Y_sh) / size_k
  foi_g  = foi_gk + (foi_Agg + foi_Ygg) / (y[gSS]*2 + y[gSA] + y[gSY])
  foi_k  = p$rate_sh_kg * (y_gA  *B_A_sh + y_gY  *B_Y_sh) / size_g +
           p$rate_sh_kk * (y[kAu]*B_A_sh + y[kYu]*B_Y_sh) / size_k
  # dy
  dy[gSu] = -y[gSu]*p$exit +y[kSu]*p$turn -form_Sx +end_Sx +p$size*p$prop_g*p$exit -y[gSu]*foi_gk +y[gSY]*p$mort
  dy[kSu] = -y[kSu]*p$exit -y[kSu]*p$turn                  +p$size*p$prop_k*p$exit -y[kSu]*foi_k
  dy[gAu] = -y[gAu]*p$exit +y[kAu]*p$turn -form_Ax +end_Ax -y[gAu]*p$prog          +y[gSu]*foi_gk +y[gAY]*p$mort
  dy[kAu] = -y[kAu]*p$exit -y[kAu]*p$turn                  -y[kAu]*p$prog          +y[kSu]*foi_k
  dy[gYu] = -y[gYu]*p$exit +y[kYu]*p$turn -form_Yx +end_Yx +y[gAu]*p$prog -y[gYu]*p$mort
  dy[kYu] = -y[kYu]*p$exit -y[kYu]*p$turn                  +y[kAu]*p$prog -y[kYu]*p$mort
  dy[gSS] = -y[gSS]*(p$exit+1/p$dur_lo) +form_SS                                                               -y[gSS]*2*foi_gk
  dy[gSA] = -y[gSA]*(p$exit+1/p$dur_lo) +form_SA -y[gSA]*p$prog                                 -y[gSA]*foi_gk +y[gSS]*2*foi_gk -foi_Agg
  dy[gSY] = -y[gSY]*(p$exit+1/p$dur_lo) +form_SY +y[gSA]*p$prog   -y[gSY]*p$mort                -y[gSY]*foi_gk -foi_Ygg
  dy[gAA] = -y[gAA]*(p$exit+1/p$dur_lo) +form_AA -y[gAA]*2*p$prog                               +y[gSA]*foi_gk +foi_Agg
  dy[gAY] = -y[gAY]*(p$exit+1/p$dur_lo) +form_AY +y[gAA]*2*p$prog -y[gAY]*p$prog -y[gAY]*p$mort +y[gSY]*foi_gk +foi_Ygg
  dy[gYY] = -y[gYY]*(p$exit+1/p$dur_lo) +form_YY -y[gYY]*2*p$mort +y[gAY]*p$prog
  return(list(dy=dy,
    foi_g  = unname(foi_g),
    foi_k  = unname(foi_k),
    prev_g = unname(prev_g),
    prev_k = unname(prev_k)))
}
# ------------------------------------------------------------------------------
# solver
solve = function(case,dt=.1,...){
  time = seq(1970,2030,dt)
  p  = par.def(...)
  dy = get(paste0('dy.',case))
  y0 = get(paste0('y0.',case))(p)
  y  = as.data.frame(ode(y0,time,dy,p))
  y[,12:17] = y[,12:17] * 2 # HACK for pair state sizes (i+1 for time)
  y$total = rowSums(y[,ylabs])
  return(y)
}
cases = c(flat='flat',epa='epa',pair='pair')
