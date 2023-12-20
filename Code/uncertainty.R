library(pracma)
library(parallel)
library(RColorBrewer)
library(tikzDevice)

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)


setwd("~/Publications/subsoilCseq/")
source("Code/BVP.R")
source("Code/GLSOM_BVP.R")
source("Code/Peclet.R")

h=0.1
d=seq(0,100, by=h)
nd=length(d)

u<-function(d, beta=0.92){
  -beta^d * log(beta)
}

ud<-u(d,beta=0.92)

k0<-1 # yr-1
efold=90 # cm
kd<-exp(-d/efold)*k0
kf=function(d,k0=1,efold=90){exp(-d/efold)*k0}

abvIn=0
x0<-(u(0)+abvIn)/k0

Ds<-runif(n, min(Kappa), max(Kappa)) #rnK
As<-runif(n, min(v), max(v)) #rnv
betas<-runif(n,0.9, 1)
ufuns<-function(bts){
  f<-function(d, beta=bts){
    -beta^d * log(beta)
  }
  return(f)
}
uds<-lapply(betas, FUN=ufuns)
kfuns<-function(efs){
  f<-function(d,k0=1,efold=efs){exp(-d/efold)*k0}
  return(f)
}
es<-runif(n, 10, 100)
kds<-lapply(es, kfuns)

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "linspace", "trisolve"))
Dunc<-clusterMap(cl, fun=GLSOM, D=Ds, MoreArgs = list(xgrid=d, a=mean(As),k=kf,f=u, boundary=c(x0,0)))
Vunc<-clusterMap(cl, fun=GLSOM, a=As, MoreArgs = list(xgrid=d, D=mean(Ds),k=kf,f=u, boundary=c(x0,0)))
Inunc<-clusterMap(cl, fun=GLSOM, f=uds, MoreArgs = list(xgrid=d, a=mean(As), D=mean(Ds),k=kf, boundary=c(x0,0)))
kunc<-clusterMap(cl, fun=GLSOM, k=kds, MoreArgs = list(xgrid=d, a=mean(As), D=mean(Ds),f=u, boundary=c(x0,0)))

uvar<-function(x){x$U}

Dxs<-sapply(Dunc, FUN=uvar)
Vxs<-sapply(Vunc, FUN=uvar)
Inxs<-sapply(Inunc, uvar)
kxs<-sapply(kunc, FUN=uvar)
Dxr<-data.frame(min=apply(Dxs, 1, min), max=apply(Dxs, 1, max))
Vxr<-data.frame(min=apply(Vxs, 1, min), max=apply(Vxs, 1, max))
Inxr<-data.frame(min=apply(Inxs, 1, min), max=apply(Inxs, 1, max))
kxr<-data.frame(min=apply(kxs, 1, min), max=apply(kxs, 1, max))

ulines<-sapply(uds, FUN=function(x){x(d)})
urange<-data.frame(min=apply(ulines, 1, min), max=apply(ulines, 1, max))
klines<-sapply(kds, FUN=function(x){x(d)})
krange<-data.frame(min=apply(klines, 1, min), max=apply(klines, 1, max))

ytm=seq(0,100, by=20)
pal<-hcl.colors(4, palette = "RdYlBu", alpha=0.5)

tikz("Figures/uncertainty.tex", standAlone = TRUE)
par(mfrow=c(2,2), mar=c(4,4,0.1, 0.1))
plot(Ds, As, ylab="$v$ (cm yr$^{-1}$)", xlab="$\\kappa$ (cm$^2$ yr$^{-1}$)", pch=20, xlim=c(0,20), ylim=c(0,10), bty="n")
points(Kappa, v, pch=20, cex=2, col=2)
legend("topleft", c("Random variates for simulation", "Values from literature"), pch=20, col=1:2, bty="n")
legend("topright", legend="a", bty="n")

plot(NA, type="l", xlim=c(0,0.2), ylim=c(-100, 0), xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", yaxt="n", lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
polygon(x=c(Vxr$max, rev(Vxr$min)), y=c(-d, rev(-d)), border=NA, col=pal[1])
polygon(x=c(Dxr$max, rev(Dxr$min)), y=c(-d, rev(-d)), border=NA, col=pal[2])
legend("bottomright", c("Uncertainty $v$", "Uncertainty $\\kappa$"), pch=15, col=pal[1:2], bty="n")
legend("topright", legend="b", bty="n")

plot(NA, type="l", xlim=c(0,1), ylim=c(-100, 0), 
     xlab=expression(paste("Input (g c", m^-2,"\\ y", r^-1, ") or decomposition (y", r^-1,")  rate")), 
     ylab="Depth (cm)", yaxt="n", lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
polygon(x=c(krange$max, rev(krange$min)), y=c(-d, rev(-d)), border = NA, col=pal[3])
polygon(x=c(urange$max, rev(urange$min)), y=c(-d, rev(-d)), border = NA, col=pal[4])
legend("bottomright", c("Uncertainty $d_e$", "Uncertainty $\\beta$"), pch=15, col=pal[3:4], bty="n")
legend("topright", legend="c", bty="n")

plot(NA, type="l", xlim=c(0,0.2), ylim=c(-100, 0), xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", yaxt="n",lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
polygon(x=c(kxr$max, rev(kxr$min)), y=c(-d, rev(-d)), border=NA, col=pal[3])
polygon(x=c(Inxr$max, rev(Inxr$min)), y=c(-d, rev(-d)), border=NA, col=pal[4])
legend("bottomright", c("Uncertainty $k(d)$", "Uncertainty $u(d)$"), pch=15, col=pal[3:4], bty="n")
legend("topright", legend="d", bty="n")
par(mfrow=c(1,1))
dev.off()

#matplot(klines, -d, type="l", lty=1, col=1, lwd=0.5)
