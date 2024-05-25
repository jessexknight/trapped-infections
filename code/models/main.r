source('models/def.r')
source('models/utils.r')

check.leak = function(){
  ym = y.melt(lapply(cases,solve,mort=0))
  plot(ym,c('total'),linetype=case)
  fig.save('leaks')
}
plot.main = function(ym,slug){
  for (var in c('prev','foi')){
    g = plot(ym,paste0(var,c('_g','_k')),linetype=case) +
      facet_grid('variable',scales='free') +
      expand_limits(y=0)
    fig.save(slug,var,w=4,h=4)
  }
}
run.fit = function(){
  # config
  p.adj    = c('rate_sh_kk','freq_lo')
  var.targ = c('prev_g','prev_k')
  y.targ   = solve('pair',dt=1)[c('time',var.targ)]
  # run fitting
  p.fit = list(
    epa  = fit('epa', y.targ,p.adj),
    flat = fit('flat',y.targ,p.adj))
}

# check.leak()
# ym.base = y.melt(lapply(cases,solve))
# plot.main(ym.base,'base')

print(run.fit())
ym.fit = y.melt(list(
  pair = solve('pair'),
  epa  = solve('epa', freq_lo=114),
  flat = solve('flat',freq_lo=37)
))
plot.main(ym.fit,'fit')
