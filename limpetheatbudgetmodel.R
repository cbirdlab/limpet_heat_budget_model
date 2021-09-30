rm(list=ls())
#from Denny and Harley 2006

#meteorological inputs#########################################################
Isw<-1   #solar irradiance (W/m^2)  ; got number from wikipedia
Ta<-25 + 273.15    #air temp in K = T in C + 273
Tw<-15 + 273.15    #water temp in K = T in C + 273
u<-2      #wind speed (m/s)
phi<-1.2 #angle of light relative to spire of limpet, 0 is directly overhead, pi/2 is from the side

#measurements from limpets######################################################
R <- 1.76/100               #radius cm /100 = m
H <- 2*R*.25                #height cm /100 = m
Acd <- pi*R^2               #Conductive area m2, area of foot, btw,opihi can control this
Al <- pi*R*(H^2 + R^2)^0.5  #Lateral area m2

Vol <- pi*R^2*H*(1/3)
Kg <- Vol*1023            #1023 kg/m^3 is density of seawater at 26C
JperC <- Kg*4010          #amount of energy to raise limpet temp by 1C if limpet is sea water with no loss, 4010 is specific heat of seawater at 26C in J/kg

if(H*sin(phi) <= R*abs(cos(phi))) { 
  Ap<-Acd*abs(cos(phi))                 #area of limpet's shell (m^2) projected in the direction at which sunlight strikes shell
}else{
  cot <- cos(phi)/sin(phi)
  ustar <- (R/H)*(H^2-R^2*cot^2)^0.5
  ah <- (2*H*sin(phi)) * ustar - ustar^3 * (H*sin(phi))/R^2 - abs(cos(phi)) * (R^2*sin(ustar/R)^-1 + ustar * (R^2-ustar^2)^0.5)
  Ap<- (Acd*abs(cos(phi)))+(ah)
}

asw<-0.68  #Short-wave absorptivity
alws<-0.97  #Long-wave absorptivity of shell
Elws<-0.97 #Long-wave emissivity of shell
Elwa<-9.2e-6*Ta^2   #clear sky   Long-wave emissivity of air, Long wave flux; Eq from Denny&Harley2006, by way of Campbell & Norman 1998 resulted in E<0.1. Mendoza et al report much higher values (up to 0.85), I altered this eq from 9.2e-7*Ta^2 to 9.2e-6*Ta^2 
#Elwa<-1            #cloudy sky   Long-wave emissivity of air, Denny&Harley 2006 say 1, but that can't be correct
#Elwb<-0.72   #emissivity of basalt

a <- 1.304 #a and b are the coefficients for the Nusselt-Reynolds number relationship (Eqn·27).
b <- 0.404  #a and b are the coefficients for the Nusselt-Reynolds number relationship (Eqn·27).

#measurements from rocks#######################################################
Kr<-3.06  #Thermal conductivity of granite rock at HMS W·m-1·K-1
rho<-2601   #Density of granite rock at HMS kg·m-3
cr<-789    #Specific heat of granite rock at HMS J·kg-1·K-1
k<-Kr/(rho*cr)   #Thermal diffusivity of granite rock at HMS m2·s-1
Vs<-0.7    #view factor, proportion of picture taken by fisheye lens perpendicular to substratum that is sky 

#measurements from air########################################################
v<- 1.25*10^-5 + 9.2e-8 * Ta   #kinematic viscosity of air
Ka<-0.00501 + 7.2e-5 * Ta      #conductivity of air
sigma<-5.6697e-8               #the Stefan-Boltzmann constant (5.6710-8·W·m-2·K-4).

#values calculated from limpet, air, and meteorological values######################

#convective heat flux
Re<-(2*u*R)/v      #Reynolds number
Nu<- a*Re^b     # (2*hc*R)/Ka    #Nusselt number
hc<- a*(Ka*u/v)^b * R^(b-1) #convective heat transfer coefficient


Wsw<-Ap*asw*Isw    #watts gained due to short wave solar irradiance
q1<-Wsw
JperC/Wsw   # this is how many secs it'll take for body temp to rise 1C without considering losses

q2<-Vs*Al*Elws*sigma*Ta^4*(Elwa-1)
q3<-4*Vs*Al*Elws*sigma*Ta^3

Acv<-Al           #area of shell in convective contact with air
q4<-hc*Acv

#jumpstart the rock temp gradient
TidalRange<-2 #tidal range in meters, the distance to rock that is same temp as ocean
deltaZ=0.01     #1 cm increments 
q5<-(Acd*Kr)/deltaZ
Tb<-Ta  #initialize body temp
Tr<-c(Tb,rep(Tw,each=(TidalRange/deltaZ)-1))  #this is the profile of temperatures in the rock between limpet and ocean
for(i in 1:10000){
  Ti<-diff(Tr)/deltaZ
  Ti2<-diff(Ti)/deltaZ
  Tr<-Tr[-(1)]
  Tr<-Tr[-(199)]
  Tr<-c(Tb,Tr+k*Ti2*30,Tw)
  T2<-Tr[2]
  Tb<-(q1+q2+q3*Ta+q4*Ta+q5*T2)/(q3+q4+q5)
}
Tr
Tb-273.15

#evaluate fluxes####################################################################
Wsw<-Ap*asw*Isw    #watts gained due to short wave solar irradiance
Wsw
Wlws<-Vs*Al*Elws*sigma*Tb^4  #watts lost due to long wave to sky
Wlws
Wlwa <-Vs*Al*alws*Elwa*sigma*Ta^4  #watts gained by long wave from sky
Wlwa
Wlw <- Wlwa - Wlws   #the net rate of long-wave energy transfer
Wlw
Wlw<-q2+q3*(Ta-Tb)
Wlw

Wcv<-q4*(Ta-Tb)   #convective heat transfer
Wcv

Wcd<-q5*(T2-Tb)   #Conductive heat flux into or out of the limpet
Wcd

# Tr<-c(Tb,rep(Tw,each=199))
# for(i in 1:20000){
#   Ti<-diff(Tr)/0.01
#   Ti2<-diff(Ti)/0.01
#   Tr<-Tr[-(1)]
#   Tr<-Tr[-(199)]
#   Tr<-c(Tb,Tr+k*Ti2*30,Tw)
# }
# Tr

