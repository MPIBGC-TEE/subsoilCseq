library(tikzDevice)

hbars<-function(y, m, sd, cl){
  arrows(x0=m-sd, y0=y, x1=m+sd, angle=90, code=3, length = 0.1, col=cl)
}

figa<-read.csv("Data/Data_figa.csv")
figb<-read.csv("Data/Data_figb.csv")
dcm<-seq(-20, -100, by=-20)
yt<-c("0-20", "20-40", "40-60", "60-80", "80-100")
age<-subset(figa,Variable=="Age")
tt<-subset(figa,Variable=="Tau")
bms<-unique(figb$Biome)

pal<-hcl.colors(2, palette = "Red-Blue")
pal2<-hcl.colors(n=length(bms), palette = "RdYlBu")

tikz(file="Figures/biomeTT.tex", standAlone = TRUE, height=7*2)
par(mar=c(4,5,1,0), mfrow=c(2,1))
plot(age$mean, dcm, type="b", ylab="Depth (cm)", xlab="Mean age and transit time (yr)", 
     xlim=c(0,5000),yaxt="n", col=pal[1], lwd=3, bty="n")
hbars(y=dcm, m=age$mean, sd=age$se, cl=pal[1])
axis(side=2,at=dcm, labels=yt, las=1)
points(tt$mean, dcm, type="b", col=pal[2], lwd=3)
hbars(y=dcm, m=tt$mean, sd=tt$se, cl=pal[2])
legend("topright", c("Age", "Transit time"), pch=1, lty=1, lwd=3, col=pal, bty="n")

plot(figb[figb$Biome==bms[1],"mean"], dcm, type="b", ylab="Depth (cm)", xlab="Mean transit time (yr)", 
     log="x", xlim=c(10,100000), xaxt="n",yaxt="n", col=pal2[1], lwd=3, bty="n")
hbars(y=dcm, m=figb[figb$Biome==bms[1],"mean"], sd=figb[figb$Biome==bms[1],"se"], cl=pal2[1])
axis(side=2,at=dcm, labels=yt, las=1)
axis(side=1, at=10^seq(1,5), labels = c("10", "100", "1000", "10000", "100000"))
for(i in 2:length(bms)){
  points(figb[figb$Biome==bms[i],"mean"], dcm, type="b", col=pal2[i], lwd=3)
  hbars(y=dcm, m=figb[figb$Biome==bms[i],"mean"], sd=figb[figb$Biome==bms[i],"se"], cl=pal2[i])
}
legend("topright", legend=bms, lty=1, pch=1, col=pal2, lwd=3, bty="n")
par(mfrow=c(1,1))
dev.off()
