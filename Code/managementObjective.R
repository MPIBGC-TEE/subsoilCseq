library(tikzDevice)
library(pracma)
library(RColorBrewer)
pal1<-brewer.pal(4,"Set1")

setwd("~/Publications/subsoilCseq/")
source("Code/BVP.R")
source("Code/GLSOM_BVP.R")

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

ref<-GLSOM(xgrid=d, D=1, a=1,k=kf,f=u, boundary=c(x0,0))
opt<-GLSOM(xgrid=d, D=1, a=100,k=kf,f=u, boundary=c(x0,0.05))

ytm=seq(0,100, by=20)

tikz(file="Figures/managementObjective.tex", standAlone = TRUE)
par(mar=c(4,4,0,0), las=1)
plot(ref$U, -d, type="l", xlim=c(0,0.1), xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", yaxt="n", col=pal1[1], lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
lines(opt$U[-((nd-10):nd)], -d[-((nd-10):nd)], col=pal1[2], lwd=3)
#abline(v=opt$U[1], lty=2)
arrows(0.01,-50, 0.055, -50, lwd=2)
legend(x=0,y=-10, legend="$\\frac{\\partial x}{\\partial d} \\ll 0$", text.col=pal1[1], cex=2, bty="n")
legend(x=0.065, y=-40, legend="$\\frac{\\partial x}{\\partial d} \\approx 0$", text.col=pal1[2], cex=2, bty="n")
legend(x=0.02, y=-50, legend="$\\frac{g(d)}{v} \\to 0$", cex=2, bty="n")
dev.off()

# Graphical abstract
pal2<-c("white", rgb(5,5,1, maxColorValue = 255))

tikz(file="GraphicalAbstract/managementObjective.tex", bg=rgb(204, 186, 150, maxColorValue = 255))
par(mar=c(4,4,0,0), las=1)
plot(ref$U, -d, type="l", xlim=c(0,0.1), xlab="Carbon stock ($x$)", xaxt="n", 
     ylab="Depth ($d$, cm)", yaxt="n",  col=pal2[1], lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
axis(side=1, labels=FALSE)
lines(opt$U[-((nd-10):nd)], -d[-((nd-10):nd)], col=pal2[2], lwd=3)
arrows(0.01,-50, 0.055, -50, lwd=2)
legend(x=0,y=-10, legend="$\\frac{\\partial x}{\\partial d} \\ll 0$", text.col=pal2[1], cex=2, bty="n")
legend(x=0.065, y=-40, legend="$\\frac{\\partial x}{\\partial d} \\approx 0$", text.col=pal2[2], cex=2, bty="n")
legend(x=0.01, y=-60, legend="$\\frac{g(d)}{v} \\to 0$", cex=2, bty="n")
legend(x=0.01, y=-50, legend="Increase subsoil C", cex=2, bty="n")
dev.off()

