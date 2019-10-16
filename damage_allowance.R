
#!/usr/local/bin/R
#
# damage_allowance.R
# Written by D.J. Rasmussen (dmr2-AT-princeton-DOT-edu)
#
# Last updated: Tue Oct 15 10:18:51 PDT 2019
# Calculate instantaneous damage allowances for the flood mitigation strategies of: 

# 1. Storm surge barrier
# 2. Levee
# 3. Coastal Retreat
# 4. Elevation...

# ...under different assumptions of 1) AIS behavior and 2) climate forcing scenario

# Input: 
#        1) R routines listed below and the R libraries: pracma, lhs
#        2) Monte Carlo SLR samples from LocalizeSL in ASCII text format (in centimeters)
#        3) GPD parameters for an extreme sea level distribution
#           Note: parameters for tide gauges around the world are provided in the included table
#        4) Property value aggregated at different elevations for constructing a damage function
#           Note: Assumes property is in billions of USD; Elevation is in meters and is relative to MHHW
#           Property value for NYC is provided in file. 
#           The depth-damage functions currently assume a 70/30 residential/commercial highrise split

# Output: a table with flood damage allowances for a given 
# damage mitigation strategy considering: 
#
#         1) the entire SLR distribution 
#         2) The 5th SLR percentile
#         3) The 95th SLR percentile
#


rm(list=ls(all=TRUE))

#project directory
root <- "/Users/dmr/Dropbox\ (Princeton)/Projects/Damage\ Allowances"
subProj <- "damage\ allowance"
setwd(paste(root,"/",subProj,sep=""))

source("routines/func.R")
source("routines/GPDsample.R")
source("routines/GPDNExceed.R")
source("routines/dmg_allwnc_routines.R")
source("routines/openSLR_noFilter.R")

library(pracma) # for 'trapz'

# Weights and thresholds we want to plot
alphaList <- c(1750,1500,1000,500,250) # 2100 max AIS contributions in millimeters
betaList <- c(0,0.25,0.5,0.75,1) # weights for AIS melt "speed"

scen = "rcp85" # climate forcing scenario

yr1 <- 2000 # Choose time frame
yr2 <- 2100


# Input files 

# Monte Carlo SLR samples
dir <- paste("data/slr_du_rasmussen19/lslr_nyc_deep_uncertainty_rasmussen19/",scen,"/",sep="")

# Flood distribution parameters for NYC (Battery Park)
scale <- 0.126937112 # median scale parameter
shape <- 0.189351584 # median shape parameter
threshold <- 0.510255768 # GPD threshold
lambda <- 2.631559987
shapescaleV <- -0.000515828 # covariance of GPD scale and shape parameters
shapeV <- 0.004965107 # variance of shape parameter
scaleV <- 0.000142102 # variance of scale parameter

# Existing Flood protection level for NYC (Battery Park)
protectLvl <- 1.06 # meters above MHHW (NOAA/NWS)

zmax <- 10.0 # max extreme water level to evaluate (meters)
dz <- 0.01 
z <- seq(0,zmax,dz) # generate some extreme water levels (relative to MHHW)


# 1-D aggregate damage function for NYC

# Read file with property information 

fil <- "data/manhattan_propval_area_10cm.csv"
df_prop <- read.csv(fil)

# fit piece-wise function to accumulated property curve
bp <- 3
z_prop <- df_prop$Elevation

fit <- lm(df_prop$prop_acc_sum[z_prop<=bp]~poly(z_prop[z_prop<=bp],2,raw=T)-1)
a1 <- fit$coefficients[1]
b1 <- fit$coefficients[2]

fit <-lm(df_prop$prop_acc_sum[z_prop>bp & z<=zmax] ~ z_prop[z_prop>bp & z_prop<=zmax])
a2 <- fit$coefficients[1]
b2 <- fit$coefficients[2]


# Get modeled marginal property function from fits
damage <- depth_dmg(z, bp, a1, b1, b2)
damage[z<protectLvl] <- 0

# Targeted annual average loss
targ_aal <- NA # billions of USD$ or 'NA' for using current AAL


# Generate Parameter Samples for Sampling Parameter Uncertainty...
GPDsamps <- GPDsample(1000, scale, shape, shapeV, scaleV, shapescaleV)

         
for (strategy in c("elevation2","levee","surge_barrier","coastal_retreat")){

# Choose compliance rates for elevation and coastal retreat
elev_compliance_rate <- 1.0 # perfect compliance
retreat_compliance_rate <- 1.0


##### END USER INPUTS #####


if (strategy=="coastal_retreat"){
  coastal_retreat = TRUE
  elevation1 = FALSE
  levee = FALSE
  surge_barrier = FALSE
  elevation2 = FALSE
}else if(strategy=="elevation1"){
  elevation1 = TRUE
  coastal_retreat = FALSE
  levee = FALSE
  surge_barrier = FALSE
  elevation2 = FALSE
}else if(strategy=="levee"){
  design_tol <- .10 # failure rate at design height
  leveeFB <- .5 # freeboard in meters
  levee = TRUE
  elevation1 = FALSE
  coastal_retreat = FALSE
  surge_barrier = FALSE
  elevation2 = FALSE
}else if(strategy=="surge_barrier"){
# Choose threshold for surge barrier closure
  zclose <- 1.0 # gate closure threshold (m above MHHW)
  design_tol <- .10 # failure rate at design height
  leveeFB <- .5 # freeboard in meters
  surge_barrier = TRUE
  elevation1 = FALSE
  levee = FALSE
  coastal_retreat = FALSE
  elevation2 = FALSE
}else if(strategy=="elevation2"){
  elevation2 = TRUE
  surge_barrier = FALSE
  elevation = FALSE
  levee = FALSE
  coastal_retreat = FALSE
  elevation1 = FALSE
}else{
  stop(paste("Invalid Strategy:",strategy))
}


if (is.na(targ_aal)){
  
  print(" **** No Flood Risk Target Found. Using Current Annual Avg. Loss... **** ")
  
  # Calculate current annual average loss
  qqq <- matrix(NaN, nrow=length(GPDsamps[,2]), ncol=length(z))
  for(iii in 1:length(GPDsamps[,2]) ){
    qqq[iii,] <- GPDNExceed(z-threshold,lambda,-threshold,GPDsamps[iii,2],GPDsamps[iii,1])
  }
  
  Nz <- apply(qqq,2,mean,na.rm=T)
  fz <- diff(exp(-apply(qqq,2,mean,na.rm=T))/dz)
  fz <- c(fz,fz[length(z)-1])
  
  targ_aal <- trapz(z,fz*damage)
}

print(paste("Using Flood Risk Target: $",round(targ_aal,3)," billion USD",sep=""))

inst <- array( NA, dim=c(length(betaList), length(alphaList), 3) ) # Instantaneous Damage Allowance


# each year 
for (yr in seq(yr1,yr2,10)){

# each beta
for (bbb in seq_along(betaList) ){
    
# each alpha
for (aaa in  seq_along(alphaList)) {
  
  if (yr == 2000){
    inst[,,] <- 0.0
    break
  }

 # Open local SLR data (assumes centimeters)
  prfx <- "LSLproj_MC_K14_K17_"
  sffx <- paste("alpha",alphaList[aaa],"mm_beta",deci2lab(betaList[bbb]),"_NEW_12_",scen,".tsv",sep="")
  fSLR <- paste(dir,prfx,sffx,sep="")
  SLR <- getRSLsamps_noFilter( fSLR )
    
  years <- SLR$years 
  indx = which(SLR$years==yr) 
  slr <- SLR$samples[,indx]*10 # cm to meters
    
  slr5 <- round( quantile(slr,probs=c(.05)), 2)
  slr95 <- round( quantile(slr,probs=c(.95)), 2)

  print(paste("The expected/5/95 SLR for ",yr," is: ",round(mean(slr),2),"m/",slr5,"m/",slr95,"m (alpha=",alphaList[aaa],"mm)",sep=""))
  
  fSLR <- array( NaN, dim=c(3, length(z)) )
  
  # DL SLR PDF
  fSLR[1,] <- ewl_wSLR_pdf(z, slr, lambda, threshold, scale, shape, shapeV, scaleV, shapescaleV)
  fSLR[2,] <- ewl_wSLR_pdf(z, slr5, lambda, threshold, scale, shape, shapeV, scaleV, shapescaleV)
  fSLR[3,] <- ewl_wSLR_pdf(z, slr95, lambda, threshold, scale, shape, shapeV, scaleV, shapescaleV)
  
if ( coastal_retreat ){
      cat("\n \n")
      print(paste("Calculating Instantaneous Coastal Retreat Damage Allowance...",sep=""))
      print(paste("Target Year:",yr,sep=" "))
      print(paste("Coastal Retreat Compliance Rate: ",retreat_compliance_rate," \\ SLR beta=",betaList[bbb],"; SLR alpha=",alphaList[aaa],"mm",sep=""))
    
      retreat_allwnc_dl <- dl_retreat(z, damage, fSLR[1,], targ_aal, retreat_compliance_rate, protectLvl) - protectLvl
      print(paste("Instantaneous Coastal Retreat Damage Allowance: ",retreat_allwnc_dl," meters [ Year: ",yr," ]",sep=""))
      inst[bbb,aaa,1] <- retreat_allwnc_dl
      inst[bbb,aaa,2] <- dl_retreat(z, damage, fSLR[2,], targ_aal, retreat_compliance_rate, protectLvl) - protectLvl
      inst[bbb,aaa,3] <- dl_retreat(z, damage, fSLR[3,], targ_aal, retreat_compliance_rate, protectLvl) - protectLvl
} # Retreat

if ( elevation1 ){
      cat("\n \n")
      print(paste("Calculating Instantaneous Elevation Damage Allowance...",sep=""))
      print(paste("Target Year:",yr,sep=" "))
      print(paste("Elevation Compliance Rate: ",elev_compliance_rate,sep=""))

      elev_allwnc_dl <- dl_elevate_all(z,bp,a1,b1,b2,fSLR[1,],targ_aal,elev_compliance_rate,protectLvl)
      print(paste("Instantaneous Elevation Damage Allowance: ",elev_allwnc_dl," meters [ Year: ",yr," ]",sep=""))
      inst[bbb,aaa,1] <- elev_allwnc_dl
      inst[bbb,aaa,2] <- dl_elevate_all(z,bp,a1,b1,b2,fSLR[2,],targ_aal,elev_compliance_rate,protectLvl)
      inst[bbb,aaa,3] <- dl_elevate_all(z,bp,a1,b1,b2,fSLR[3,],targ_aal,elev_compliance_rate,protectLvl)
} # Elevation
  
if ( elevation2 ){
      cat("\n \n")
      print(paste("Calculating Instantaneous Elevation Damage Allowance (below A) ...",sep=""))
      print(paste("Target Year:",yr,sep=" "))
      print(paste("Elevation Compliance Rate: ",elev_compliance_rate,sep=""))
      
      elev_allwnc_dl <- dl_elevate_belowA(z,bp,a1,a2,b1,b2,fSLR[1,],targ_aal,elev_compliance_rate,protectLvl) - protectLvl
      print(paste("Instantaneous Elevation Damage Allowance: ",elev_allwnc_dl," meters [ Year: ",yr," ]",sep=""))
      inst[bbb,aaa,1] <- elev_allwnc_dl
      inst[bbb,aaa,2] <- dl_elevate_belowA(z,bp,a1,a2,b1,b2,fSLR[2,],targ_aal,elev_compliance_rate,protectLvl) - protectLvl
      inst[bbb,aaa,3] <- dl_elevate_belowA(z,bp,a1,a2,b1,b2,fSLR[3,],targ_aal,elev_compliance_rate,protectLvl) - protectLvl
} # Elevation

# Levee 
if ( levee ){
      cat("\n \n")
      print(paste("Calculating Instantaneous Levee Damage Allowance...",sep=""))
      print(paste("Target Year:",yr,sep=" "))
      print(paste("Design Failure Tolerance:",design_tol,"// Levee Freeboard:",leveeFB,"meters",sep=" "))

      levee_allwnc_dl <- dl_levee(z, damage, fSLR[1,], targ_aal, design_tol, leveeFB) - protectLvl
      print(paste("Instantaneous Levee Damage Allowance: ",levee_allwnc_dl," meters [ Year: ",yr," ]",sep=""))
      inst[bbb,aaa,1] <- levee_allwnc_dl
      inst[bbb,aaa,2] <- dl_levee(z, damage, fSLR[2,], targ_aal, design_tol, leveeFB) - protectLvl
      inst[bbb,aaa,3] <- dl_levee(z, damage, fSLR[3,], targ_aal, design_tol, leveeFB) - protectLvl
} # Levee

if ( surge_barrier ){
      cat("\n \n")
      print(paste("Calculating Instantaneous Surge Barrier Damage Allowance...",sep=""))
      print(paste("Target Year:",yr,sep=" "))
      print(paste("Design Failure Tolerance:",design_tol,"// Barrier Freeboard:",leveeFB,"meters",sep=" "))
    
      barrier_allwnc_dl <- dl_barrier(z, damage, fSLR[1,], targ_aal, design_tol, leveeFB, zclose) - protectLvl
      print(paste("Instantaneous Surge Barrier Damage Allowance: ",barrier_allwnc_dl," meters [ Year: ",yr," ]",sep=""))
      inst[bbb,aaa,1] <- barrier_allwnc_dl
      inst[bbb,aaa,2] <- dl_barrier(z, damage, fSLR[2,], targ_aal, design_tol, leveeFB, zclose) - protectLvl
      inst[bbb,aaa,3] <- dl_barrier(z, damage, fSLR[3,], targ_aal, design_tol, leveeFB, zclose) - protectLvl
} # Surge Barrier
  
} # each alpha
} # each beta
 

# Write Results to Disk
row.names(inst) <- betaList
colnames(inst) <- alphaList

# Expected SLR (full distribution)
outf <- paste("saved_results/dl_allwnc_expected_slr_du_",strategy,"_",scen,"_",yr,".tsv",sep="")
write.table(inst[,,1],file=outf, sep="\t", quote = FALSE,col.names=NA)

# 5th percentile SLR
outf <- paste("saved_results/dl_allwnc_5th_slr_du_",strategy,"_",scen,"_",yr,".tsv",sep="")
write.table(inst[,,2],file=outf, sep="\t", quote = FALSE,col.names=NA)

# 95th  percentile SLR
outf <- paste("saved_results/dl_allwnc_95th_slr_du_",strategy,"_",scen,"_",yr,".tsv",sep="")
write.table(inst[,,3],file=outf, sep="\t", quote = FALSE,col.names=NA)

} # each year
} # each strategy
