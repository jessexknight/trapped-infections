library('ggplot2')
library('reshape2')
# infer project root folder
proj.root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
root.path = function(...,create=FALSE){
  # define & maybe create a path from root
  path = file.path(proj.root,...)
  if (create){ dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE) }
  return(path)
}
y.melt = function(ylist){
  # melt a list of model outputs & add case column
  ym = do.call(rbind,lapply(names(ylist),function(case){
    yi = cbind(melt(ylist[[case]],id='time'),case=case)
  }))
}
plot = function(ym,vars,...){
  # plot melted model outputs
  yp = ym[ym$variable %in% vars,]
  g = ggplot(yp,aes(x=time,y=value,...)) +
    geom_line() +
    theme_light()
}
fig.save = function(...,w=4,h=3){
  # save fig in proj.root/fig/
  fname = paste(...,'pdf',sep='.')
  fpath = root.path('fig',fname,create=TRUE)
  ggsave(fpath,w=w,h=h)
}
