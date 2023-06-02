#' A General Linear Second Order Model for representing depth profiles at equilibrium
#' 
#' This function implements a general linear second order model for representing vertical transport 
#' at steady-state. It uses a second order differencing method to solve for the equilibrium state.
#' given some boundary conditions. 
#' @param x A vector with the spatial coordinates of the vertical profile
#' @param h The size of the intervals in x. A scalar value
#' @param D Diffusion coefficient. A scalar value
#' @param a Advection rate. A scalar value
#' @param c Reaction rate. A scalar value
#' @param f A vector with the vertical input rate
#' @param boundary A vector of length 2 with the upper and lower boundary conditions. These are Neumann boundary conditions with the input rate on the top layer as first element, and the bottom output rate as second element

GLSOM<-function(x, h, D, a, c, f, boundary){
  
  m<-length(x)
  alpha<-boundary[1]
  beta<-boundary[2]
  
  # Original formula from Leveque assuming Dirichlet boundary conditions at both ends
  F<-matrix(f,nrow=m, ncol=1)
  F[1]<-F[1]-((D/(h^2)) - (a/(2*h)))*alpha
  F[m]<-F[m]-((D/(h^2)) + (a/(2*h)))*beta
  
  # Build matrix A as a linear combination of D, a, and c.
  AD<-diag(-2/(h^2),nrow=m)
  AD[row(AD)-1 == col(AD)]<- 1/(h^2)
  AD[row(AD) == col(AD)-1]<- 1/(h^2)

  # AD<-diag(0,nrow=m)
  # AD[1,1:5]<-c(2.9167e+04,  -8.6667e+04,   9.5000e+04,  -4.6667e+04,   9.1667e+03)
  # AD[2,1:5]<-c(9.1667e+03,  -1.6667e+04,   5.0000e+03,   3.3333e+03,  -8.3333e+02)
  # AD[3,1:5]<-c(-8.3333e+02,   1.3333e+04,  -2.5000e+04,   1.3333e+04,  -8.3333e+02)
  # for(i in 1:(m-5)){
  #   AD[(3+i),(1+i):(5+i)]<-c(-8.3333e+02,   1.3333e+04,  -2.5000e+04,   1.3333e+04,  -8.3333e+02)
  # }
  # AD[(m-1),(m-6):(m-1)]<-c(-8.3333e+03,   5.0833e+04,  -1.3000e+05,   1.7833e+05,  -1.2833e+05,   3.7500e+04)
  # AD[m,(m-5):m]<-c(-8.3333e+03,   5.0833e+04,  -1.3000e+05,   1.7833e+05,  -1.2833e+05,   3.7500e+04)
  
  # B<-diag(0,nrow=m)
  # B[row(B)-1 == col(B)]<- -1/(2*h)
  # B[row(B) == col(B)-1]<- 1/(2*h)
  # B[1,1:2]<-c(-100,100) # Fornberg coefficients
  

  # # Fourth order approximation for h = 0.01
  # B<-diag(0,nrow=m)
  # B[1,1:5]<-c(-208.333,   400.000,  -300.000,   133.333,   -25.000)
  # B[2,1:5]<-c(-25.0000,   -83.3333,   150.0000,   -50.0000,     8.3333)
  # B[3,1:5]<-c(8.3333e+00,  -6.6667e+01,  0,   6.6667e+01,  -8.3333e+00)
  # B[4,1:5]<-c(-8.3333,    50.0000,  -150.0000,    83.3333,    25.0000)
  # for(i in 1:(m-5)){
  #   B[(4+i),(1+i):(5+i)]<-c(-8.3333,    50.0000,  -150.0000,    83.3333,    25.0000)
  # }
  # B[m,(m-5):m]<-c(-20.000,   125.000,  -333.333,   500.000,  -500.000,   228.333)

  B<-diag(0,nrow=m)
  B[1,1:5]<-c(1/4, -4/3,  3,    -4,    25/12)/h
  B[2,1:5]<-c(1/4, -4/3,  3,    -4,    25/12)/h
  B[3,1:5]<-c(1/4, -4/3,  3,    -4,    25/12)/h
  B[4,1:5]<-c(1/4, -4/3,  3,    -4,    25/12)/h
  for(i in 1:(m-5)){
    B[(4+i),(1+i):(5+i)]<-c(1/4, -4/3,  3,    -4,    25/12)/h
  }
  B[m,(m-4):m]<-c(1/4, -4/3,  3,    -4,    25/12)/h
  
  A<- D*AD + a*B + diag(c, nrow=m)
  
#   # build only tridiagonal elements and use function trisolve
#   d<- rep((-2/(h^2))+c, m)
# #  d[1]<-d[1]+(-100) #Fornberg weighs second order
#   ud<- rep((1/(2*h))+(1/(h^2)), m-1)
# #  ud[1]<-ud[1]+50 #Fornberg weighs second order
#   ld<- rep((-1/(2*h))+(1/(h^2)), m-1)
#   
#   
#   A<-diag(d)
#   A[row(A)-1 == col(A)]<-ld
#   A[row(A) == col(A)-1]<-ud
#   
#  U<-trisolve(a=d,b=ld,d=ud,rhs=F)
  
  U<- solve(A)%*%F
  Out<-list(U, A, F) # pack as list
  names(Out)<-c("U", "A", "F")
  return(Out)
}


