
vElzein1995<-c(0.13, 0.34, 0.48, 0.6, 0.42) # mm yr-1
kappaElzein1995<-c(5.15, 16.58, 5.29, 0.94, 1.48) # cm2 yr-1

vBruun2007<-0.0081 # cm yr-1
kappaBruun2007<-0.71 # cm2 yr-1 

vBraakhekke2011<-0.002 # m yr-1
kappaBraakhekke2011<- 0.4*0.3/(1000*2) *10000 # cm2 yr-1. Obtained as B*lm/(rho*2). The last term changes units from m2 yr-1 to cm2 yr-1.
                          # rho is bulk density, and here a value of 1000 kg cm-3 is assumed.

# Values from page 106, Braakhekke's thesis. First value corresponds to Lobos, second value is for Hainich, Mode B (most likely parameter combination)
vBraakhekke2013<-c(0.0651, 0.00137) # m yr-1
kappaBraakhekke2013<-c(0.00943*2.53/(1400*2) *10000, # Lobos. cm2 yr-1. Obtained as B*lm/(rho*2)
                       0.233*0.583/(1000*2) * 10000)  # Hainich, Mode B. We assume a bulk density of 1000 kg cm-3

v<-c(vElzein1995/10, # divide by 10 cm to chage to cm yr-1
     vBruun2007,
     vBraakhekke2011*100, # multyply by 100 cm to get cm yr-1
     vBraakhekke2013*100 # multyply by 100 cm to get cm yr-1
)
Kappa<-c(kappaElzein1995,
         kappaBruun2007,
         kappaBraakhekke2011,
         kappaBraakhekke2013
         )

Pe<-v/Kappa
round(Pe, 3)

hist(Pe)
boxplot(Pe, log="y")
abline(h=1)

plot(v, Kappa, xlim=c(0,16), ylim=c(0,16), pch=19, xlab="Advection velocity", ylab="Diffusion coefficient", bty="n")
abline(0,1, lty=2)
segments(-1,5,5,5, col=2)
segments(5,5, 5, -1, col=2)
points(c(0.25, 0.5, 0.5), c(0.5, 0.5, 0.25), pch=19, col=2)
