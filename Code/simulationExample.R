library(SoilR)
library(expm)
library(RColorBrewer)
library(tikzDevice)
library(knitr)
library(pracma)
# library(extrafont)
# font_import()
# loadfonts()
pal1<-brewer.pal(4,"Set1")
pal2<-brewer.pal(4,"Dark2")

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

plot(ud, -d, type="l")
plot(kd, -d, type="l")

Ds<-c(0.1,1,5,15)
As<-c(0.1,1,5,10)

Druns<-mapply(FUN=GLSOM,D=Ds, MoreArgs = list(xgrid=d, a=As[2],k=kf,f=u, boundary=c(x0,0)))
Vruns<-mapply(FUN=GLSOM,a=As, MoreArgs = list(xgrid=d, D=Ds[2],k=kf,f=u, boundary=c(x0,0)))

ytm=seq(0,100, by=20)

#pdf("Figures/simulationsKappaV.pdf")
tikz(file="Figures/simulationsKappaV.tex", standAlone = TRUE)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), cex=1.1)
plot(Druns[[1]], -d, type="l", xlim=c(0,0.1), xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", yaxt="n", col=pal1[1], lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
lines(Druns[[4]], -d, col=pal1[2], lwd=3)
lines(Druns[[7]], -d, col=pal1[3], lwd=3)
lines(Druns[[10]], -d, col=pal1[4], lwd=3)
legend("bottomright", legend=Ds, title=expression(paste("$\\kappa$", "\\ (c", m^2, "\\ y", r^-1, ")")), lty=1, lwd=3, col=pal1, bty="n")

plot(diff(Druns[[1]])[-((nd-100):nd)], -d[-((nd-100):nd)], xlim=c(-0.001,0.001), ylim=c(-100,0), 
     xaxt="n", yaxt="n", type="l", xlab=expression(paste("First derivative of concentration (x", 10^-3,")")), 
     ylab=" ", col=pal1[1], lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
axis(side=1, at=seq(-1e-03, 1e-03, 5e-04), labels=seq(-1,1, by=0.5))
lines(diff(Druns[[4]]), -d[-1], col=pal1[2], lwd=3)
lines(diff(Druns[[7]]), -d[-1], col=pal1[3], lwd=3)
lines(diff(Druns[[10]]), -d[-1], col=pal1[4], lwd=3)
abline(v=0, lty=2, lwd=1)

plot(Vruns[[1]][-nd], -d[-nd], type="l", xlim=c(0,0.1), xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", yaxt="n", col=pal2[1], lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
lines(Vruns[[4]][-nd], -d[-nd], col=pal2[2], lwd=3)
lines(Vruns[[7]][-nd], -d[-nd], col=pal2[3], lwd=3)
lines(Vruns[[10]][-((nd-20):nd)], -d[-((nd-20):nd)], col=pal2[4], lwd=3)
legend("bottomright", legend=As, title=expression(paste(italic(v),"\\ (cm y",r^-1,")")), lty=1, col=pal2, lwd=3, bty="n")

plot(diff(Vruns[[1]])[-nd], -d[-nd], type="l", xlim=c(-0.001,0.001), xaxt="n",
     xlab=expression(paste("First derivative of concentration (x", 10^-3,")")), 
     ylab=" ", yaxt="n", col=pal2[1], lwd=3, bty="n")
axis(side=2,at=-ytm, labels=ytm)
axis(side=1, at=seq(-1e-03, 1e-03, 5e-04), labels=seq(-1,1, by=0.5))
lines(diff(Vruns[[4]])[-((nd-100):nd)], -d[-((nd-100):nd)], col=pal2[2], lwd=3)
lines(diff(Vruns[[7]])[-((nd-100):nd)], -d[-((nd-100):nd)], col=pal2[3], lwd=3)
lines(diff(Vruns[[10]])[-((nd-100):nd)], -d[-((nd-100):nd)], col=pal2[4], lwd=3)
abline(v=0, lty=2, lwd=1)
par(mfrow=c(1,1))
dev.off()
####################################

ks<-c(1,0.1, 0.01)

kdruns<-list()
for(i in 1:length(ks)){
  kdruns[[i]]<-GLSOM(xgrid=d,D=0.01,a=0.01,k=function(d){kf(d,k0=ks[i])},f=u, boundary=c(u(0)/ks[i],0))$U
}

betas<-c(0.92, 0.95, 0.98)

udruns<-list()
for(i in 1:length(betas)){
  udruns[[i]]<-GLSOM(xgrid=d,D=0.01,a=0.01,k=function(d){kf(d,k0=ks[2])},f=function(d){u(d,beta=betas[i])}, boundary=c(u(0,betas[i])/ks[2],0))$U
}

#pdf("Figures/rootDecomp.pdf")
tikz(file="Figures/rootDecomp.tex", standAlone = TRUE)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), lwd=2, cex=1.1)
plot(kdruns[[3]], -d, type="l", yaxt="n", xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", col=pal1[1], lwd=3, bty="n")
axis(side=2, at=-ytm, labels=ytm)
lines(kdruns[[2]], -d, col=pal1[2], lwd=3)
lines(kdruns[[1]], -d, col=pal1[3], lwd=3)
legend("bottomright", legend=rev(ks),title=expression(paste(italic(k[0]), "\\ (y", r^-1, ")")), lty=1, lwd=3, col=pal1, bty="n")

plot(diff(kdruns[[3]])[-((nd-20):nd)], -d[-((nd-20):nd)], type="l", yaxt="n", xlab="First derivative of concentration", ylab=" ", col=pal1[1], lwd=3, bty="n")
axis(side=2, at=-ytm, labels=ytm)
lines(diff(kdruns[[2]])[-((nd-20):nd)], -d[-((nd-20):nd)], col=pal1[2], lwd=3)
lines(diff(kdruns[[1]])[-((nd-20):nd)], -d[-((nd-20):nd)], col=pal1[3], lwd=3)

plot(udruns[[1]], -d, type="l", yaxt="n", xlab=expression(paste("C concentration (g c",m^-3,")")), 
     ylab="Depth (cm)", col=pal2[1], lwd=3, bty="n")
axis(side=2, at=-ytm, labels=ytm)
lines(udruns[[2]][-((nd-20):nd)], -d[-((nd-20):nd)], col=pal2[2], lwd=3)
lines(udruns[[3]][-((nd-20):nd)], -d[-((nd-20):nd)], col=pal2[3], lwd=3)
legend("bottomright", legend=betas, title="$\\beta$", lty=1, col=pal2, lwd=3, bty="n")

plot(diff(udruns[[1]])[-((nd-30):nd)], -d[-((nd-30):nd)], type="l", yaxt="n", xlab="First derivative of concentration", ylab=" ", col=pal2[1], lwd=3, bty="n")
axis(side=2, at=-ytm, labels=ytm)
lines(diff(udruns[[2]])[-((nd-30):nd)], -d[-((nd-30):nd)], col=pal2[2], lwd=3)
lines(diff(udruns[[3]])[-((nd-30):nd)], -d[-((nd-30):nd)], col=pal2[3], lwd=3)
par(mfrow=c(1,1))
dev.off()


##############
# Transport-decomposition runs

ThDf<-GLSOM(xgrid=d,D=1,a=5,k=function(d){kf(d,k0=ks[1])},f=function(d){u(d,betas[2])}, boundary=c(u(0,betas[2])/kf(0,k0=ks[1]),0))
ThDs<-GLSOM(xgrid=d,D=1,a=5,k=function(d){kf(d,k0=ks[2])},f=function(d){u(d,betas[2])}, boundary=c(u(0,betas[2])/kf(0,k0=ks[2]),0))
TlDf<-GLSOM(xgrid=d,D=1,a=0.1,k=function(d){kf(d,k0=ks[1])},f=function(d){u(d,betas[2])}, boundary=c(u(0,betas[2])/kf(0,k0=ks[1]),0))
TlDs<-GLSOM(xgrid=d,D=1,a=0.1,k=function(d){kf(d,k0=ks[2])},f=function(d){u(d,betas[2])}, boundary=c(u(0,betas[2])/kf(0,k0=ks[2]),0))

udm<-u(d,betas[2])[c(-1,-nd)]

tau=seq(0,10, by=0.1)
TT1<-transitTime(A=(h^-2)*ThDf$A, u=udm, q=0.5, a=tau)
TT2<-transitTime(A=(h^-2)*ThDs$A, u=udm, q=0.5, a=tau)
TT3<-transitTime(A=(h^-2)*TlDf$A, u=udm, q=0.5, a=tau)
TT4<-transitTime(A=(h^-2)*TlDs$A, u=udm, q=0.5, a=tau)

save(TT1,file="Code/TT1.RData")
save(TT2,file="Code/TT2.RData")
save(TT3,file="Code/TT3.RData")
save(TT4,file="Code/TT4.RData")

SA1<-systemAge(A=(h^-2)*ThDf$A, u=udm, a=tau[1], q=0.5) # Only interested in mean and median ages, so tau = 0.
SA2<-systemAge(A=(h^-2)*ThDs$A, u=udm, a=tau[1], q=0.5)
SA3<-systemAge(A=(h^-2)*TlDf$A, u=udm, a=tau[1], q=0.5)
SA4<-systemAge(A=(h^-2)*TlDs$A, u=udm, a=tau[1], q=0.5)

plot(SA1$meanPoolAge, -dm, type="l", xlim=c(0,50))
lines(SA2$meanPoolAge, -dm, col=2)
lines(SA3$meanPoolAge, -dm, col=3)
lines(SA4$meanPoolAge, -dm, col=4)

#pdf("Figures/transitTimes.pdf")
tikz(file="Figures/transitTimes.tex", standAlone = TRUE, height=7*sqrt(2))
par(mar=c(4,4,0,0), mfrow=c(2,1))
plot(tau, TT1$transitTimeDensity, type="l", log="y", ylim=c(0.001, 1), xlim=c(0,10), yaxt="n", xlab="Transit time (yr)", 
     ylab="Probability density", col=pal2[1], lwd=3, bty="n")
axis(side=2, at=10^(seq(-4, 1, by=1)), labels = c("0.0001", "0.001", "0.01", "0.1", "0", "1"))
lines(tau, TT2$transitTimeDensity, col=pal2[2], lwd=3)
lines(tau, TT3$transitTimeDensity, col=pal2[3], lwd=3)
lines(tau, TT4$transitTimeDensity, col=pal2[4], lwd=3)
legend("topright", "a", bty="n")

plot(ThDf$U[-((nd-30):nd)], -d[-((nd-30):nd)], type="l", xlim=c(0,0.8), lwd=3,
     xlab=expression(paste("C concentration (g c",m^-3,")")), yaxt="n",ylab="Depth (cm)", col=pal2[1], bty="n")
axis(side=2, at=-ytm, labels=ytm)
lines(ThDs$U[-((nd-50):nd)], -d[-((nd-50):nd)], col=pal2[2], lwd=3)
lines(TlDf$U[-((nd-30):nd)], -d[-((nd-30):nd)], col=pal2[3], lwd=3)
lines(TlDs$U[-((nd-30):nd)], -d[-((nd-30):nd)], col=pal2[4], lwd=3)
legend("bottomright", legend=c("Transport fast, Decomposition fast",
                               "Transport fast, Decomposition slow",
                               "Transport slow, Decomposition fast",
                               "Transport slow, Decomposition slow"), lty=1, lwd=3, col=pal2, bty="n")
legend("topright", "b", bty="n")
par(mfrow=c(1,1))
dev.off()

TDlist<-list(ThDf,ThDs,TlDf,TlDs)

Mfunc<-function(t, B, u){
  expm::expm(t*B)%*%(u/sum(u))
}

Ts<-c(1,5,10,50)
MThDf<-sapply(Ts,FUN=function(x){Mfunc(t=x, B=(h^-2)*ThDf$A, u=udm)})
MThDs<-sapply(Ts,FUN=function(x){Mfunc(t=x, B=(h^-2)*ThDs$A, u=udm)})
MTlDf<-sapply(Ts,FUN=function(x){Mfunc(t=x, B=(h^-2)*TlDf$A, u=udm)})
MTlDs<-sapply(Ts,FUN=function(x){Mfunc(t=x, B=(h^-2)*TlDs$A, u=udm)})

MTs<-array(0,dim=c(dim(MThDf),4))
MTs[,,1]<-MThDf; MTs[,,2]<-MThDs
MTs[,,3]<-MTlDf; MTs[,,4]<-MTlDs

dm<-d[c(-1,-nd)]
m<-length(dm)

#pdf("Figures/Mt.pdf")
tikz(file="Figures/Mt.tex", standAlone = TRUE)
par(mfrow=c(2,2), mar=c(4,4,1,0.5))
matplot(x=MTs[-((m-20):m),1,], y=-dm[-((m-20):m)], type="l", xlim=c(0,max(MTs[,1,])), lty=1, lwd=3, yaxt="n", 
        xlab=expression(paste("Proportion remaining after 1 year (c",m^-1,")")), ylab="Depth (cm)", col=pal2, bty="n")
axis(side=2,at=-ytm, labels=ytm)
matplot(x=MTs[-((m-20):m),2,], y=-dm[-((m-20):m)], type="l", xlim=c(0,max(MTs[,1,])), lty=1, lwd=3, yaxt="n", 
        xlab=expression(paste("Proportion remaining after 5 years (c",m^-1,")")), ylab="", col=pal2, bty="n")
axis(side=2,at=-ytm, labels=ytm)
matplot(x=MTs[-((m-20):m),3,], y=-dm[-((m-20):m)], type="l", xlim=c(0,max(MTs[,1,])), lty=1, lwd=3, yaxt="n", 
        xlab=expression(paste("Proportion remaining after 10 years (c",m^-1,")")), ylab="Depth (cm)", col=pal2, bty="n")
axis(side=2,at=-ytm, labels=ytm)
matplot(x=MTs[-((m-20):m),4,], y=-dm[-((m-20):m)], type="l", xlim=c(0,max(MTs[,1,])), lty=1, lwd=3, yaxt="n", 
        xlab=expression(paste("Proportion remaining after 50 years (c",m^-1,")")), ylab="", col=pal2, bty="n")
axis(side=2,at=-ytm, labels=ytm)
legend("bottomright", legend=c("Transport fast, Decomposition fast",
                              "Transport fast, Decomposition slow",
                              "Transport slow, Decomposition fast",
                              "Transport slow, Decomposition slow"), lty=1, lwd=3, col=pal2, bty="n")
par(mfrow=c(1,1))
dev.off()

sapply(TDlist, FUN=function(x){sum(x$U)})

TTlist<-list(TT1,TT2, TT3, TT4)

MTT<-sapply(TTlist, function(x) x$meanTransitTime)
QTT<-sapply(TTlist, function(x) x$quantiles)

xss<-sapply(TDlist, FUN=function(x){sum(solve(x$A)%*%x$F)})

SAlist=list(SA1,SA2,SA3,SA4)

MSA<-sapply(SAlist, function(x) x$meanSystemAge)
QSA<-sapply(SAlist, function(x) x$quantiles)


# Table
expTable<-data.frame(Parameter=c("kappa", "v","k0", "beta", "M(1)", "M(10)", "M(50)", "Mean transit time", "Median transit time", "CS(infty)"),
           TfDf=round(c(1,5,ks[1], betas[2], sum(MThDf[,1]), sum(MThDf[,3]), sum(MThDf[,4]), MTT[1], QTT[1], xss[1]),3),
           TfDs=round(c(1,5,ks[2], betas[2], sum(MThDs[,1]), sum(MThDs[,3]), sum(MThDs[,4]), MTT[2], QTT[2], xss[2]),3),
           TsDf=round(c(1,0.1,ks[1], betas[2], sum(MTlDf[,1]), sum(MTlDf[,3]), sum(MTlDf[,4]), MTT[3], QTT[3], xss[3]),3),
           TsDs=round(c(1,0.1,ks[2], betas[2], sum(MTlDs[,1]), sum(MTlDs[,3]), sum(MTlDs[,4]), MTT[4], QTT[4], xss[4]),3)
           )
kable(expTable, format="latex")

###
# Green functions obtained from Lanczos (1997, Dover) Eq. 5.16.1 (page 248)
G1<-TriSolve(B=(h^-2)*(ThDf$A), u=-u(d, betas[2])[c(-1, -nd)])
G2<-TriSolve(B=(h^-2)*(ThDs$A), u=-u(d, betas[2])[c(-1, -nd)])
G3<-TriSolve(B=(h^-2)*(TlDf$A), u=-u(d, betas[2])[c(-1, -nd)])
G4<-TriSolve(B=(h^-2)*(TlDs$A), u=-u(d, betas[2])[c(-1, -nd)])

Gs=data.frame(G1,G2,G3, G4)
matplot(d[c(-1,-nd)], Gs, type="l", lty=1, col=pal2)
legend("topright", legend=c("Transport fast, Decomposition fast",
                            "Transport fast, Decomposition slow",
                            "Transport slow, Decomposition fast",
                            "Transport slow, Decomposition slow"), lty=1, lwd=3, col=pal2, bty="n")
