
#############################################################################################################
## R Codes to calculate 
## (1) mean sound pressure levels for each octave and signals associated with flow sound/turbulence and surfacing bubbles
## (2) bubble size distributions from ambient underwater sound spectra and 
## (3) bubble residence and equilibration times
##
## Code is written by Marcus Klaus
## Swedish University of Agricultural Science
## marcus.klaus@posteo.net
## last update: 2021-06-24
##
## if you find any mistakes, I would appreciate if you let me know about them!
##
## This code is used in the following manuscript and should be cited as follows:
## Marcus Klaus, Thierry Labasque, Gianluca Botter, Nicola Durighetto, Jakob Schelker (submitted): 
## Unravelling the contribution of turbulence and bubbles to air-water gas exchange in running waters
#############################################################################################################

###########################
## load libraries
###########################

library(seewave)
library(emdbook)
library(tuneR)
library(Hmisc)

###########################
## load support functions
###########################

######################################################################################
## compute kinematic viscosity of water for given temperature
######################################################################################

## code by Jordan S. Read, from GDopp package
## available at https://rdrr.io/github/USGS-R/GDopp/src/R/get_kin_viscosity.R

get_kin_viscosity <- function(temperature=20) {
  # from Mays 2005, Water Resources Engineering
  tempTable <- seq(0,100,by=5)
  # table in m2/s E-6
  visTable <- c(1.792,1.519,1.308,1.141,1.007,0.897,
                0.804,0.727,0.661,0.605,0.556,0.513,0.477,0.444,
                0.415,0.39,0.367,0.347,0.328,0.311,0.296)
  v <- data.frame(approx(tempTable,visTable,xout = temperature))[2]
  v <- v*1e-6
  return(v$y)
}

######################################################################################
# compute freshwater density (Kg/m3)   
######################################################################################

# function from
# Chen, C. T., and F. J. Millero. 1977. The use and misuse of pure water 
# PVT properties for lake waters. Nature 266: 707–708. doi:10.1038/266707a0

get_waterdensity <- function(Temp=20) {
density = 0.999835 + (6.7914e-5*Temp) - (9.0894e-6*Temp^2) + (1.0171e-7*Temp^3) - ( 1.2846e-9   *  Temp^4) + ( 1.1592e-11  *  Temp^5) - ( 5.0125e-14  *  Temp^6) + S * (8.221e-4 - ( 3.87e-6  *  Temp) + ( 4.99e-8  * Temp^2 ))
ph=density*1000
return(ph)
}

########################################################################################################################################
## (1) load function to compute mean sound pressure levels of octaves and signals related to flow sound/turbulence and bubbles
########################################################################################################################################
# mean sound pressure levels for each octave as in 
# Klaus, M., E. Geibrink, E. R. Hotchkiss, and J. Karlsson. 2019. Listening to air–water gas exchange in running waters. Limnol. Oceanogr. Methods 17: 395–414. 
# doi:10.1002/lom3.10321
#
# mean sound pressure levels for flow sound / turbulence as in 
# Geay, T., P. Belleudy, C. Gervaise, H. Habersack, J. Aigner, A. Kreisler, H. Seitz, and J. B. Laronne. 2017.
# Passive acoustic monitoring of bed load discharge in a large gravel bed river.
# J. Geophys. Res. Earth Surf. 122: 528–545. doi:10.1002/2016JF004112
########################################################################################################################################
# INPUTS:
# SigWave = signal wave spectrum (with bubble sound)
# WaveDir = file name / directory of sound wave
# pref = reference sound pressure level (10^-6 Pa for underwater sound)
#
# OUTPUTS:
# table with mean+-sd sound pressure levels ("SPLm" and "SPLsd") for octaves from 16 Hz to 16 KHz, 
# and for intervals associated to flow sound / turbulence and surfacing bubbles (prmst and prmsb)
########################################################################################################################################

dbf <- function(SigWave,WaveDir,pref=10^-6){

db <- data.frame(matrix(0,1,25))
colnames(db) <- c("file","SPLm16","SPLm31.5","SPLm63","SPLm125","SPLm250","SPLm500","SPLm1k","SPLm2k","SPLm4k","SPLm8k","SPLm16k","SPLsd16","SPLsd31.5","SPLsd63","SPLsd125","SPLsd250","SPLsd500","SPLsd1k","SPLsd2k","SPLsd4k","SPLsd8k","SPLsd16k","Prmst(Pa)","PrmsB (Pa)")

### means
db[,1] <- WaveDir
db[,2] <- mean(SigWave[SigWave[,1]>=0.011 & SigWave[,1]<0.022,2])
db[,3] <- mean(SigWave[SigWave[,1]>=0.022 & SigWave[,1]<0.044,2])
db[,4] <- mean(SigWave[SigWave[,1]>=0.044 & SigWave[,1]<0.088,2])
db[,5] <- mean(SigWave[SigWave[,1]>=0.088 & SigWave[,1]<0.177,2])
db[,6] <- mean(SigWave[SigWave[,1]>=0.177 & SigWave[,1]<0.355,2])
db[,7] <- mean(SigWave[SigWave[,1]>=0.355 & SigWave[,1]<0.710,2])
db[,8] <- mean(SigWave[SigWave[,1]>=0.710 & SigWave[,1]<1.42,2])
db[,9] <- mean(SigWave[SigWave[,1]>=1.42 & SigWave[,1]<2.84,2])
db[,10] <- mean(SigWave[SigWave[,1]>=2.84 & SigWave[,1]<5.68,2])
db[,11] <- mean(SigWave[SigWave[,1]>=5.68 & SigWave[,1]<11.36,2])
db[,12] <- mean(SigWave[SigWave[,1]>=11.36 & SigWave[,1]<22.72,2])

# root mean square pressure (power of signal) related to turbulence (Pa) between 10 and 100 Hz (Geay et al. 2017)
# here, I exclude the interval 44-57 Hz because it was dominated by background noise (distinct peak at 50 Hz)
intt <- which((SigWave[,1]>0.01 & SigWave[,1]<0.044) | (SigWave[,1]>0.057 & SigWave[,1]<0.1))
db[,24] <- 2*sum(SigWave[intt,2]*c(diff(SigWave[intt,1])[1],diff(SigWave[intt,1])))*pref

# root mean square pressure (power of signal) related to surface bursts of bubbles (Pa), take "id" interval from bubble size model above
intsb <- which((SigWave[,1]> 5 & SigWave[,1]<10) )
db[,25] <- 2*sum(SigWave[intsb,2]*c(diff(SigWave[intsb,1])[1],diff(SigWave[intsb,1])))*pref

### standard deviations
db[,13] <- sd(SigWave[SigWave[,1]>=0.011 & SigWave[,1]<0.022,2])
db[,14] <- sd(SigWave[SigWave[,1]>=0.022 & SigWave[,1]<0.044,2])
db[,15] <- sd(SigWave[SigWave[,1]>=0.044 & SigWave[,1]<0.088,2])
db[,16] <- sd(SigWave[SigWave[,1]>=0.088 & SigWave[,1]<0.177,2])
db[,17] <- sd(SigWave[SigWave[,1]>=0.177 & SigWave[,1]<0.355,2])
db[,18] <- sd(SigWave[SigWave[,1]>=0.355 & SigWave[,1]<0.710,2])
db[,19] <- sd(SigWave[SigWave[,1]>=0.710 & SigWave[,1]<1.42,2])
db[,20] <- sd(SigWave[SigWave[,1]>=1.42 & SigWave[,1]<2.84,2])
db[,21] <- sd(SigWave[SigWave[,1]>=2.84 & SigWave[,1]<5.68,2])
db[,22] <- sd(SigWave[SigWave[,1]>=5.68 & SigWave[,1]<11.36,2])
db[,23] <- sd(SigWave[SigWave[,1]>=11.36 & SigWave[,1]<22.72,2])

return(db)
}


#########################################################################################
## (2) function to convert underwater sound spectra to bubble size spectra 
## method following Loewen and Melville (1991) J. Acoust. Soc. Am. 90(4)2075-2080
## 
## Loewen, M. R., and W. K. Melville. 1991. A model of the sound generated by breaking waves. 
## J. Acoust. Soc. Am. 90: 2075–2080. doi:https://doi.org/10.1121/1.401634
##
## Please note the assumptions behind this conversion: 
## (1) the sound spectrum is only due to bubbles oscillating 
## at their natural frequencies in a free field
## (2)single bubble's contribution
## to the spectrum is only in a narrow frequency band about f
## (3) bubble pressure pulses are damped following a sinosoidal law
## (4) bubble pressure pulses of neighboring bubbles do not affect each other  
## 
#########################################################################################
#
# INPUT VARIABLES 
# SigWaveDir = directory to signal wave file (with bubble sound)
# RefWaveDir = directory to reference wave file (without bubble sound)
# channel = hydrophone channel (here, select "L" for left and "R" for right
# Temp = Water temperature (°C)
# S = total dissolved salt (parts per thousand)
# P0 = pressure of surrounding liquid (Pa)
# d = distance between hydrophone and the free water surface (m)
# L = twice the distance between the bubble center and the free water surface (m); here assumed to be the water column depth
# R = distance between hydrophone receiver and bubble (m) 
# Qb = volumetric bubble flux (L/min)
# ovlp = overlap between two successive fft analysis windows (in %) (seewave package)
# wl = length of the window for the fft analysis (seewave package)
# Diff = subtract reference wave with background noise from signal wave? (TRUE/FALSE)
# plot = plot bubble size and sound spectra? (TRUE/FALSE)
# Filter = filter spectrum by running median to remove spikes? (TRUE/FALSE)
# nfilt = filter window size
#
# CONSTANTS
# ga = ratio of the specific heats of the bubble gas (adiabatic compression), for air it is 1.4
# ep = amplitude of the bubble surface oscillations divided by the equilibrium bubble radius
#
# CALCULATED VARIABLES
# ph = Water density (Kg/m3)
# c = speed of sound (m/s) 
# P = mean square of measured sound spectrum level (dB re 1 µPa)
# fg = frequencies associated to P (1/s)
# f = selected frequency (1/s)
# de = bubble damping coefficient (dimensionless)
# k = acoustic wave number (m)
# ag = bubble radius associated to f
# p2 = mean square signal level for a single bubble resonating at frequency f (Pa^2)
#
# OUTPUT 
# list with 5 elements:
# [[1]] power spectral density of sound pressure levels and bubble size distributions of singal wave (with bubbles)
# [[2]] power spectral density of sound pressure levels of reference wave (without bubbles)
# [[3]] summary of input and output data
# [[4]] sound pressure levels of signal wave as means for 10 octaves and for frequency bands associated with flow sound/turbulence and surfacing bubbles (see code (1))
# [[5]] sound pressure levels of reference wave channel as means for 10 octaves and for frequency bands associated with flow sound/turbulence and surfacing bubbles (see code (1))
#
# SELECTED OUTPUT VARIABLES (for details, see code below)
# a = bubble radius (m)
# N = relative number of bubbles oscillating at frequency f
# Nb = number flux of bubbles (min-1)
# Qbh = number flux of bubbles (s-1) and per nozzle of the linear aeration system used in Klaus et al (submitted)
# amean = average bubble radius weighted by number of bubbles
# asd = standard deviation of bubble radius weighted by number of bubbles
##################################################################### 

bubbleSize <- function(
SigWaveDir=SigWaveDir,
RefWaveDir=RefWaveDir,
channel="L",
Temp=10,
S=0.11,
P0=101325,
ep=0.015,
d=0.035,
L = 0.1,
R=0.05,
Qb=250,
ga=1.4,
ovlp=50,
wl = 8192*2,
id = id,
Diff=TRUE,
plot=TRUE,
Filter=FALSE,
nfilt=15,
pref=pref){


# compute freshwater density (Kg/m3) using Chen and Millero (1977), Nature 266:707-708   
ph=get_waterdensity(Temp)

## compute speed of sound (m/s) (Medwin 1975; The Journal of the Acoustical Society of America 58, 1318)
c = 1449.2 + 4.6*Temp-0.055*Temp^2+0.00029*Temp^3 + (1.34-0.010*Temp)*(S-35)+0.016*d

## load signal sound wave
SigWave <- lapply(SigWaveDir, readWave)
SigWave_R <- mono(SigWave[[1]], which = "right")
SigWave_L <- mono(SigWave[[1]], which = "left")

# set start and stop time, exclude first and last second
start=2
Sigstop=duration(SigWave_L)-2

# calculate sound pressure spectrum (power spectral density)
if(channel=="L"){
SigWave <- meanspec(SigWave_L,wl = wl, ovlp = ovlp,PSD=TRUE,PMF=FALSE,norm=FALSE,from=start,to=Sigstop,plot=FALSE)
}
if(channel=="R"){
SigWave <- meanspec(SigWave_R,wl = wl, ovlp = ovlp,PSD=TRUE,PMF=FALSE,norm=FALSE,from=start,to=Sigstop,plot=FALSE)
}
if(Filter==TRUE){
SigWave[,2] <- runmed(SigWave[,2],nfilt)
}

## load reference sound wave
RefWave <- lapply(RefWaveDir, readWave)
RefWave_R <- mono(RefWave[[1]], which = "right")
RefWave_L <- mono(RefWave[[1]], which = "left")
Refstop=duration(RefWave_L)-2

# calculate sound pressure spectrum (power spectral density)
if(channel=="L"){
RefWave <- meanspec(RefWave_L,wl = wl, ovlp = ovlp,PSD=TRUE,PMF=FALSE,norm=FALSE,from=start,to=Refstop,plot=FALSE)
}
if(channel=="R"){
RefWave <- meanspec(RefWave_R,wl = wl, ovlp = ovlp,PSD=TRUE,PMF=FALSE,norm=FALSE,from=start,to=Refstop,plot=FALSE)
}
if(Filter==TRUE){
RefWave[,2] <- runmed(RefWave[,2],nfilt)
}

# convert frequency from kHz to Hz
f=SigWave[,1]*1000
# convert f to radian frequency (rad/s)
w=2*pi*f

# bubble damping coefficient (energy loss due to thermal conduction, work against viscous forces + energy lost as acoustic wave)
de = (4.4*10^-4*f^0.5)/(1+f/(2.5^10^5))  # Crowther (1988), cited in Loewen and Melville (1991)

# acoustic wave number (1/m)
k = 2*pi*f/c

## bubble radius associated to f
a = 1/(w)*(3*ga*P0/ph)^(1/2)

#################################
## compute bubble size spectrum

# get measured sound power spectrum density (square of spectrum) (Pa2/Hz)
# if Diff==TRUE, subtract reference (background noise) from signal
if(Diff==TRUE){
P <- SigWave[,2]-RefWave[,2]
}else{
P <- SigWave[,2]
}

## mean square signal level of single bubble resonating at frequency f (Pa2)
p2 =  (3*ga*P0/ph)^3*(ph*ep*d*L/(R^2*c))^2* ( (2/(de*w*(de^2+4))+(de^2+2)/((k*R)^2*(de*w*(de^2+4))))-2/(k*R*w*(de^2+4) ) ) 

#  number of bubbles oscillating at frequency f
N = (P/p2)

# compute N only if the following assumptions are fulfilled:
# - k*L<1, this removes extremely small bubbles 
# - k*a should be <<1
# note that the maximum bubble radius found in other studies was 7.4 mm (Loewen and Melville 1991)
# bubbles are likely smaller here because my spectra decrease rapidly beyond around 6 mm

id = which(k*L<1 & k*a<0.05) 

## number flux of bubbles (min-1) = (m3/min/m3)
Nb <- matrix(NA,length(N),1)
Nb[id] <- Qb/1000*N[id]/(sum(N[id]*(4/3*pi*a[id]^3),na.rm=T))

## number flux of bubbles (s-1) per noozle of the linear aeration system 
# (for a 35 m long linear aeration system with nozzles in 2 rows and 0.013 m distance inbetween)
Qbh = sum(Nb,na.rm=T)/60/(2*35/0.013)

# number average bubble size and measure of spread (following 
# Pandit, A. B., J. Varley, R. B. Thorpe, and J. F. Davidson. 1992. 
# Measurement of Bubble Size Distribution: An acoustic technique. 
# Chem. Eng. Sci. 47: 1079–1089. doi:https://doi.org/10.1016/0009-2509(92)80233-3

amean = sum(Nb[id]*a[id],na.rm=TRUE)/sum(Nb[id],na.rm=TRUE)
Pmean = sum(Nb[id]*P[id],na.rm=TRUE)/sum(Nb[id],na.rm=TRUE)

# standard deviation of bubble size
asd = sqrt(wtd.var(a[id], Nb[id]))
Psd = sqrt(wtd.var(P[id], Nb[id]))

##################################################
### plot bubble size spectra and sound spectra
##################################################

if(plot==TRUE){

plot(a*1000,Nb,type="l",ylab=expression(paste("# Bubbles [ ",flume^-1, min^-1,"]")),xlab="Bubble radius [mm]",
ylim=c(0,15000),xlim=c(0,8),col=8)
abline(v=amean*1000)
lines(a[id]*1000,Nb[id],col=1,lwd=2)


local({
dev.set (2)
dev.print (device=jpeg, file=paste0(gsub(".wav","",SigWaveDir),"_Chan",channel,"_BubbleSize.jpg"), width=par("din")[1]*100, res=100);
})
dev.off()

plot(SigWave[,1],SigWave[,2],log="xy",type="l",ylim=c(1,max(SigWave[,2],na.rm=T)),col=8,xlab="f (kHz)",ylab="PSD (Pa2/Hz)")
lines(RefWave[,1],RefWave[,2],col="orange")
abline(v=c(0.016,0.031,0.063,0.125),lty=2)

lines(SigWave[id,],lwd=2,col=1)
lines(RefWave[id,],lwd=2,col=2)

local({
dev.set (2)
dev.print (device=jpeg, file=paste0(gsub(".wav","",SigWaveDir),"_Chan",channel,"_SpectraSigRef.jpg"), width=par("din")[1]*100, res=100);
})
dev.off()
}

#########################################################################################################
## extract sound pressure for each octave and at frequencies related to flow/turbulence and surfacing bubbles
#########################################################################################################

dbSig <- dbf(SigWave,SigWaveDir,pref=pref)
dbRef <- dbf(RefWave,RefWaveDir,pref=pref)

# add variable that is TRUE if conditions for bubble size estimation are fulfilled
BubbleTest <- k*L<1 & k*a<0.05

####################################
### make output files
####################################

SigSpectrum <- data.frame(SigWaveDir,f,SigWave[,2],a,BubbleTest,Nb,N,P)
RefSpectrum <- data.frame(RefWaveDir,f,RefWave[,2])

colnames(SigSpectrum) <- c("FileID","f Frequency (Hz)","SigWave PSD (Pa2/Hz)","a Radius (m)","Valid bubbles","Nb #bubbles (min-1)",
"N rel #bubbles","P PSD (Pa2/Hz)")

colnames(RefSpectrum) <- c("FileID","f Frequency (Hz)","P PSD (Pa2/Hz)")

metadata <- data.frame(SigWaveDir,Temp,S,P0,d,L,R,ep,Qb,Qbh,amean,asd,Pmean,Psd,Filter,nfilt,Diff)

return(list(SigSpectrum,RefSpectrum,metadata,dbSig,dbRef))
}



##############################################################################
##############################################################################
### (3) function to estimate bubble residence time and equilibration time
### of a bubble released at depth z and rising in cross flow with velocity u
##############################################################################
##############################################################################
## this code is largely following equations in 
# Woolf, D. K. 1993. Bubbles and the air - sea transfer velocity of gases Bubbles and the Air-Sea Transfer Velocity of Gases. 
# Atmosphere-Ocean 31: 517–540. doi:10.1080/07055900.1993.9649484
#
## INPUTS
# a = bubble radius (m)
# g = gravity accelaration (m/s2)
# u = flow velocity (m/s)
# Temp = water temperature (°C)
# z = water depth (m)
# Schmidt = Schmidt number of the gas
# Ostwald = Ostwald solubility of the gas in water
# gas = gas, e.g. "He", "CO2" (this text will be appended to output headers)
#
## OUTPUTS
# Tstar = dimensionless lifetime of gas in bubble
# Tg = Bubble equilibration time (s)
# T = Bubble residence time (s)
# j = gas transfer velocity of single bubble (m/s) (outputs from various models)
# Re = Reynolds number
##############################################################################

getTstar <- function(a=0.001,g=9.81,Temp=10,u=0.1,z=0.1,Schmidt=600,Ostwald=1,gas="He"){

# kinematic viscosity
kinvis = get_kin_viscosity(Temp)

# surface tension of the air-water interface (N/m) fitted to data by Jasper (1972), T is water temperature (10-25°C)
ss=(-0.148*Temp + 75.84)/1000

# water density (Kg/m3)
ph=get_waterdensity(Temp)

## vertical bubble rise (slip) velocity (m/s) of bubble in cross flow with radius > 0.65 mm; Zhang et al. (2014), eq. 11,  from Clift et al. (1978) cited therein
us=sqrt(2.14*ss/(ph*a)+0.505*g*a)

### Bubble rise velocity in cross flow (Zhang et al. 2014)(m/s)
# assuming that bubbles accelate each other, reasonable for a > 3.5 mm
# bubble induced rise velocity (m/s) (Zhang et al. 2013)
ubw = (0.94*a*2*1000-0.29)/100
# here assume that ub and ubw are additive (true additional velocity is likely slightly less than ubw)
ub = sqrt(us^2+u^2) + ubw

## distance the bubble travels in water (m)
Db = z/sin(atan(us/u))

## time the bubble is underwater (Bubble lifetime (s))
T = Db/ub

## Reynolds number [] (Memery and Merlivat 1985)
Re = 2*ub*a/kinvis

D <- NULL 
j <- matrix(0,length(a),length(gas)) 
jRe <- matrix(0,length(a),length(gas)) 
jmax <- matrix(0,length(a),length(gas)) 
jMem <- matrix(0,length(a),length(gas)) 
Tg <- matrix(0,length(a),length(gas))
Tstar <- matrix(0,length(a),length(gas))

for(i in 1:length(gas)){
# Diffusivity (m2/s)
D[i] = kinvis/Schmidt[i]

## gas transfer velocity of single bubble (m2/s*m/s /m)^1/2 = (m/s)
# The hydro- dynamical conditions of the water surrounding the bubble 
# change when the surface is assumed to be liquid (clean bubble) or solid (dirty bubble)
# T for large bubbles < 1 s, probably too short for a bubble to become dirty (Detwiler and Blanchard, 1978) 
## assume a clean bubble

# there are different models for j (see below). we here use the model recommended by Woolf 1993:

## for Re >= 10
jRe[,i] = ((1-2.89/sqrt(Re))*2*D[i]*ub/(pi*a))^(1/2)
## for Re < 10
jRe[,i][which(Re<10)] = 0.23*(D[i][which(Re<10)]*ub/a[which(Re<10)])^1/2

# alternative model by Woolf 1993, eq 43, leading to higher estimates
jmax[,i] = (1.75 * 10^-8*a^-3 + 1.84*a^-0.25)*D[i]^(1/2)

# Memery and Merlivat 1985; this model predicts much higher j than the other models!
# Memery, L., and L. Merlivat. 1985. Modelling of gas flux through bubbles at the air-water interface.
# Tellus 37B: 272–285. doi:https://doi.org/10.3402/tellusb.v37i4-5.15030
jMem[,i] <- 8*sqrt(pi*D[i]*ub/(2*a))

# use here the model recommended by Woolf 1993, but also report results from other models
j[,i]=jRe[,i]

## Bubble equilibration time (s)
Tg[,i] = a/(3*j[,i]*Ostwald[i]) 

## dimensionless lifetime of gas in bubble
Tstar[,i] = T/Tg[,i]
}

colnames(Tstar) <- paste0("Tstar_",gas)
colnames(Tg) <- paste0("Tg_",gas)
colnames(j) <- paste0("j_",gas)
colnames(jRe) <- paste0("jRe_",gas)
colnames(jmax) <- paste0("jmax_",gas)
colnames(jMem) <- paste0("jMem_",gas)

return(cbind(Tstar,T,Tg,j,jRe,jmax,jMem,Re))

}

## end of code





