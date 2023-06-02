library(pracma)
source(GLSOM_Dirilecht)
# solve y'' + y = x from Braun (1993, Springer) page 152, example 1
# with analytical solution: y(x)=cos(x)+sin(x)+x
# and y(0)=1, y(10)=8.61
h<-0.01
x<-seq(0,10,by=h)

y<-function(x) cos(x)+sin(x)+x

yt<-GLSOM(x,h,D=1,a=0,c=1,f=x,boundary = c(1,y(tail(x,1))))
yt<-BVP(f=0,g=-1,h=function(x){x},x=range(x),y=c(1,y(tail(x,1))),n=length(x)-2)
plot(x,yt$U)
lines(x,y(x),col=2, lwd=2)



##########
# Solve y'' + y' + y = x^2 from Braun (1993) page 158, example 1
# analytical solution: y(x) = -2x + x^2

x<-seq(0,10,by=h)
y2<-function(x) -2*x + x^2
y2p<-function(x) -2+2*x

yt2<-GLSOM(x,h,D=1, a=1, c=1, f=x^2, c(y2(0), y2(tail(x,1))))
yt2<-BVP(f=-1,g=-1,h=function(x){x^2},x=range(x),y=c(y2(0),y2(tail(x,1))),n=length(x)-2)
plot(x,yt2$U)
lines(x,y2(x),col=2, lwd=2)

ybvp<-bvp(f=-1, g=-1, h=function(x){x^2}, x=range(x), y=c(y2(x[1]), y2(tail(x,1))), n=length(x))
lines(as.data.frame(ybvp), col=4)

ysh<-shooting(f=function(t,y1,y2){-y2-y1+t^2},t0=x[1],tfinal=tail(x,1), y0=y2(x[1]), 
              h=function(u,v){y2(u)}, a=-2, b=-2.0)
lines(ysh$t,ysh$y[,1],col=3)
