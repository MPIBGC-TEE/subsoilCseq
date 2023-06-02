h=0.1
d=seq(0,10,by=h)
dmax=length(d)
k=0.01

x0=1
x10=0


A=diag(-2,nrow=length(d))
A[row(A)-1 == col(A)]<-1
A[row(A) == col(A)-1]<-1
B=(1/(h^2))*A

f=-exp(-d)/k
f[1]<-f[1]-(x0/(h^2))
f[dmax]<-f[dmax]-(x10/(h^2))

U=solve(B)%*%f

plot(as.numeric(U), -d)

plot(f, -d)
