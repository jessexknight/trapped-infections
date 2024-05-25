source('models/def.r')
source('models/utils.r')

check.leak = function(){
  ym = y.melt(lapply(cases,solve,mort=0))
  plot(ym,c('total'),linetype=case)
  fig.save('leaks')
}
base.case = function(){
  ym = y.melt(lapply(cases,solve))
  for (var in c('prev','foi')){
    g = plot(ym,paste0(var,c('_g','_k')),linetype=case) +
      facet_grid('variable',scales='free') +
      expand_limits(y=0)
    fig.save('base',var,w=4,h=4)
  }
}

# check.leak()
# base.case()
