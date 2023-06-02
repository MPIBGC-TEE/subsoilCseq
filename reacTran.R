library(ReacTran)

parms <- c(F0 = 1, v = 1, k = 0.1, dx = 1)

advModel <- function(t, C, parms) {
  
  with (as.list(parms), {
    
    Tran <- tran.1D(C = C, D = 0, flux.up = F0, v = v, dx = dx)
    Consumption <-  k*C
    dC   <- Tran$dC - Consumption
    
    return (list(dC = dC, Consumption= Consumption,
                 flux.up = Tran$flux.up, flux.down = Tran$flux.down))
  })
  
}

out <- steady.1D(func = advModel, y = runif(25), parms = parms,
                 nspec = 1, positive = TRUE)

plot (out, xlab = "x", ylab = "Conc", main = "advection")

ts<-seq(0,10,by=0.5)

sol=ode.1D(y=rep(0,10),times=ts,func=advModel, parms=parms,nspec = 1)

###
library(SoilR)
t=seq(from=0,to=10,by=1/12)
Litter=data.frame(year=c(1:10),Litter=rnorm(n=10,mean=10,sd=2))
years=seq(from=1,to=10,by=1/365)
TempData=data.frame(years,Temp=15+sin(2*pi*years)+
                      rnorm(n=length(years),mean=0,sd=1))
TempEffect=fT.Q10(Temp=TempData[,2],Q10=1.4)
TempEffects=data.frame(years,TempEffect)
C0=c(C01=10,C02=40,C03=30) #Initial conditions for each pool

ks=c(k1=1/2, k2=1/5, k3=1/8) #Decomposition rates for each pool
a21=ks[1]*0.3 #Transfer coefficient from pool 1 to 2
a12=ks[2]*0.1 #Transfer coefficient from pool 2 to 1
a32=ks[2]*0.2 #Transfer coefficient from pool 2 to 3
a23=ks[3]*0.1 #Transfer coefficient from pool 3 to 2
MyModel=ThreepFeedbackModel(years[1:100],ks,a21,a12,a32,a23,C0,
                            In=Litter,xi=TempEffects)


getC(MyModel)
SoilR:::getInFluxes(MyModel)@map(1)
SoilR:::getOutputFluxes(MyModel)

rso=SoilR:::getRightHandSideOfODE(MyModel)
x=rso(y=c(1,1,1), t=2)

###############################################
river.model <- function (t = 0, OC, pars = NULL) {
  
  tran <- tran.volume.1D(C = OC, F.up = F.OC, F.lat = F.lat,
                         Disp = Disp, flow = flow.up, flow.lat = flow.lat, 
                         V = Volume, full.output = TRUE) 
  
  reac <- - k*OC
  return(list(dCdt = tran$dC + reac, Flow = tran$flow))
}

nbox          <- 500     # number of grid cells
lengthEstuary <- 100000  # length of estuary [m]
BoxLength     <- lengthEstuary/nbox # [m]
Distance      <- seq(BoxLength/2, by = BoxLength, len =nbox) # [m]
Int.Distance  <- seq(0, by = BoxLength, len = (nbox+1))      # [m]

# Cross sectional area: sigmoid function of estuarine distance [m2]
CrossArea <- 4000 + 72000 * Distance^5 /(Distance^5+50000^5)

# Volume of boxes                          (m3)
Volume  <- CrossArea*BoxLength

# Transport coefficients
Disp    <- 1000   # m3/s, bulk dispersion coefficient
flow.up  <- 180    # m3/s, main river upstream inflow
flow.lat.0  <- 180    # m3/s, side river inflow

F.OC    <- 180               # input organic carbon [mol s-1]
F.lat.0 <- 180              # lateral input organic carbon [mol s-1]

k       <- 10/(365*24*3600)  # decay constant organic carbon [s-1]


#scenario 1: without lateral input
F.lat    <- rep(0, length.out = nbox)
flow.lat <- rep(0, length.out = nbox)

Conc1 <- steady.1D(runif(nbox), fun = river.model, nspec = 1, name = "OC")   

#scenario 2: with lateral input
F.lat <- F.lat.0 * dnorm(x =Distance/lengthEstuary,
                         mean = Distance[nbox/2]/lengthEstuary, 
                         sd = 1/20, log = FALSE)/nbox 
flow.lat <- flow.lat.0 * dnorm(x = Distance/lengthEstuary,
                               mean = Distance[nbox/2]/lengthEstuary, 
                               sd = 1/20, log = FALSE)/nbox 

Conc2 <- steady.1D(runif(nbox), fun = river.model, nspec = 1, name = "OC")   

plot(Conc1, Conc2, grid = Distance/1000, which = "OC", 
     mfrow = c(2, 1), lwd = 2, xlab = "distance [km]", 
     main = "Organic carbon decay in the estuary",
     ylab = "OC Concentration [mM]")
plot(Conc1, Conc2, grid = Int.Distance/1000, which = "Flow", 
     mfrow = NULL, lwd = 2, xlab = "distance [km]", 
     main = "Longitudinal change in the water flow rate",
     ylab = "Flow rate [m3 s-1]")  

legend ("topright", lty = 1:2, col = 1:2, lwd = 2,
        c("baseline", "+ side river input"))
par(mfrow=c(1,1))
####################################################
## Elzein & Balesdent 1995 Model

EB95<-function(t=0, OC, parms=NULL){
  
  tran <- tran.volume.1D(C = OC, F.up = F.OC, F.lat = F.lat,
                         Disp = Disp, flow = v, flow.lat = flow.lat, 
                         V = Volume, full.output = TRUE) 
  
  poolModel<-ThreepSeriesModel(t,ks,a21,a32,C0=rep(0,3),In=F.OC)
  rhs <-SoilR:::getRightHandSideOfODE(poolModel)
  reac <- sum(rhs(c(F.lat.0, 0, 0)))
#  reac <- - k*OC
  return(list(dCdt = tran$dC + reac, Flow = tran$flow))
  
}

ks<-c(k1=0.268, k2=0.0031, k3=1e-05) # yr-1
a21=0.05*ks[1] # yr-1
a32=0.04*ks[2] # yr-1
Disp=1.48*10e-04 # cm2 yr-1 x 10-4 to convert to m2/yr
v<-0.42*0.001 # from mm yr-1 to m yr-1

nbox          <- 100     # number of grid cells
lengthProfile <- 1  # depth of profile [m]
BoxLength     <- lengthProfile/nbox # 1 m
Distance      <- seq(BoxLength/2, by = BoxLength, len =nbox) # [m]
Int.Distance  <- seq(0, by = BoxLength, len = (nbox+1))      # [m]

# Volume of boxes                          (m3)
Volume  <- rep(0.0001,length(Distance)) # boxes of 1 cm3


fd0<-0.48 # kg m-2 yr-1 , dividing by 1 m gives kg m-3 yr-1
fd1<-3.1 # m-1
fd<-fd0*exp(-fd1*Distance/BoxLength)
F.OC<-0.17 # fs

F.lat<-fd
flow.lat<-1


soilModel <- steady.1D(fd, fun = EB95, nspec = 1, name = "OC")   

plot(soilModel)
