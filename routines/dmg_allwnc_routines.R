# Damage allowance functions
source("/Users/dmr/Dropbox\ (Princeton)/Projects/Damage\ Allowances/deep\ uncertainty/routines/GPDENExceed.R")
source("/Users/dmr/Dropbox\ (Princeton)/Projects/Damage\ Allowances/deep\ uncertainty/routines/GPDNExceed.R")
source("/Users/dmr/Dropbox\ (Princeton)/Projects/Damage\ Allowances/deep\ uncertainty/routines/GPDsample.R")
library(pracma) # for 'trapz'

ewl_wSLR_pdf <- function(z, slr, lambda, threshold, scale, shape, shapeV, scaleV, shapescaleV){
  
  dz <- diff(z)[1]
  
  # Generate Parameter Samples for Sampling Parameter Uncertainty...
  GPDsamps <- GPDsample(1000, scale, shape, shapeV, scaleV, shapescaleV)
  
  if ( length(slr) > 1 ){
    zminusSLR <- mapply(function(x) x - slr, z) - threshold
    
    # Get number of expected events/yr
    maxinterpN <- 1000
    ans <- GPDENExceed(zminusSLR, lambda, -threshold, GPDsamps[,2], GPDsamps[,1], maxinterpN)
    result <- apply( matrix(exp(ans), nrow=length(slr), ncol=1001), 2, mean, na.rm=T)
    pdf <- diff( exp(-result)/dz )
    
  } else {
    qqq <- matrix(NaN, nrow=length(GPDsamps[,2]), ncol=length(z))
    for(iii in 1:length(GPDsamps[,2]) ){
      qqq[iii,] <- GPDNExceed(z-threshold-slr,lambda,-threshold,GPDsamps[iii,2],GPDsamps[iii,1])
    }  
    # Convert number of expected events/yr to probability density function
    pdf <- diff(exp(-apply(qqq,2,mean,na.rm=T))/dz)
  }

  pdf <- c(pdf,pdf[length(z)-1])
  
  return( pdf )
}

current_aal <- function(z, damage, lambda, threshold, scale, shape, shapeV, scaleV, shapescaleV){
  
  # Generate Parameter Samples for Sampling Parameter Uncertainty...
  GPDsamps <- GPDsample(1000, scale, shape, shapeV, scaleV, shapescaleV)
  
  qqq <- matrix(NaN, nrow=length(GPDsamps[,2]), ncol=length(z))
  for(iii in 1:length(GPDsamps[,2]) ){
    qqq[iii,] <- GPDNExceed(z-threshold,lambda,-threshold,GPDsamps[iii,2],GPDsamps[iii,1])
  }
  
  fz <- diff(exp(-apply(qqq,2,mean,na.rm=T))/dz)
  fz <- c(fz,fz[length(z)-1])
  
  return( trapz(z,fz*damage) )
}

eff_pdf <- function(slr_pdf1, slr_pdf2, beta){
  
  eff_slr_pdf <- beta*sort(slr_pdf1) + (1-beta)*sort(slr_pdf2)
  return( eff_slr_pdf )
}


pfail <- function(z, prtct_hght, design_tol, leveeFB){
  
  # Estimate probability of levee failure
     b <- log(1/design_tol)/leveeFB
     a <- exp(-b*prtct_hght + log(design_tol))
    
     return( pmin(a*exp(b*z), 1) )
}

## depricated ##

dmg_atan <- function(s,protectLvl,a,b,c){
  
  dam1 <- (s+a*atan(s/b))*c
  
  # No damage below protection height
  dam1[s<=protectLvl] <- 0
  
  dam1[dam1<0] <- 0
  
  return(dam1)
} # dmg_atan

## depricated ##

dl_retreat <- function(z,damage,fSLR,target,alpha,protectLvl){
   
   A <- z
   dz <- diff(z)[1]
   
   aal_retreat <- NA
   belowA <- NA
   aboveA <- NA
  
  for (i in 1:length(A)){
     belowA[i] <- (1-alpha)*dz*sum(fSLR[1:i]*damage[1:i])
     aboveA[i] <- dz*sum(fSLR[(i+1):length(z)]*(damage[(i+1):length(z)]-alpha*damage[i]))
     aal_retreat[i] <- aboveA[i] + belowA[i]
  }

   if (max(aal_retreat,na.rm = T) < target){
     print("Cannot find an allowance for coastal retreat. The AAL target may be too large...")
     print(paste("target=",target,sep=""))
     return( NA )
   }
   
   tol <- 5e-2
   diff <- abs(aal_retreat-target)
   if ( min(diff,na.rm = T) > tol ){
     print("Didn't find a damage stabilizing protection height for Coastal Retreat!")
     print(paste("Try increasing the compliance rate...     compliance_rate=",alpha,sep=""))
     return( NA )
   }else{
     return( max(A[which.min(diff)],protectLvl) )
   }
}


dl_elevate_all <- function(z,bp,a1,b1,b2,fSLR,target,alpha,protectLvl){
  
   A <- z
   dz <- diff(z)[1]
  
   z1 <- sapply( A, function(x) z-alpha*x)
   
   # Solve for A using method of bisection
   a <- 0
   b <- 5
   
   tol <- .005
   n <- 1
   nmax <- 15
   while (n <= nmax){
     c <- (a + b)/2
     
     damage_c <- depth_dmg(z-c, bp, a1, b1, b2) 
     damage_c[z<protectLvl] <- 0
     
     damage_a <- depth_dmg(z-a, bp, a1, b1, b2) 
     damage_a[z<protectLvl] <- 0
     
     fc <- trapz(z,fSLR*damage_c) - target
     fa <- trapz(z,fSLR*damage_a) - target
     if (fc == 0 || (b-a)/2 < tol){
       
       return(c)
     }
     n <- n + 1
     
     if (sign(fc) == sign(fa)){ # new interval
       a <- c
     }else{
       b <- c
     }
   }

     print("Cannot find an allowance for elevation. The AAL target may be too large...")
     print(paste("target=",target,sep=""))
     return( NA )
}


dl_elevate_belowA <- function(z,bp,a1,a2,b1,b2,fSLR,target,alpha,protectLvl){
  
  A <- z
  dz <- diff(z)[1]
  
  z1 <- sapply( A, function(x) z-alpha*x)
  
  # Solve for A using method of bisection
  a <- 0
  b <- 5
  
  tol <- .005
  n <- 1
  nmax <- 15
  while (n <= nmax){
    c <- (a + b)/2
    
    damage_c <- depth_dmg_elev(z, c, bp, a1, b1, a2, b2)
    damage_c[z<protectLvl] <- 0
    
    damage_a <- depth_dmg_elev(z, a, bp, a1, b1, a2, b2)
    damage_a[z<protectLvl] <- 0
    
    fc <- trapz(z,fSLR*damage_c) - target
    fa <- trapz(z,fSLR*damage_a) - target
    if (fc == 0 || (b-a)/2 < tol){
      
      return(c)
    }
    n <- n + 1
    
    if (sign(fc) == sign(fa)){ # new interval
      a <- c
    }else{
      b <- c
    }
  }
  
  print("Cannot find an allowance for elevation. The AAL target may be too large...")
  print(paste("target=",target,sep=""))
  return( NA )
}

dl_levee <- function(z, damage, fSLR, target, design_tol, leveeFB){
   A <- z
   
   dz <- diff(z)[1]
   aal_levee <- NA
   for (i in 1:length(A)){ # calculate AAL under different levee heights (including free board)
     pf <- pfail(z, A[i], design_tol, leveeFB)
     aal_levee[i] <- dz*sum( fSLR[i:length(z)]*pf[i:length(z)]*damage[i:length(z)] ) 
   }

   tol <- 1e-2
   diff <- abs(aal_levee-target)
   if ( min(diff) > tol ){
     print("Could not find a damage stabilizing protection height for Levee!")
     print(paste("AAL target=",target," max AAL levee=",max(aal_levee),sep=""))
     return( NA )
   }else{
     return( A[which.min(diff)] + leveeFB )
   }
}


dl_barrier <- function(z, damage, fSLR, target, design_tol, barrierFB, zclose){
  
   dz <- diff(z)[1]
   iclose <- which.min(abs(zclose-z))
   
   # target protection height (A) must be > zclose
   A <- z[(iclose+1):length(z)]
   barrier_open <- dz*sum(fSLR[1:iclose]*damage[1:iclose])
   
   aal_barrier <- NA
   for (i in 1:length(A)){
     
     pf <- pfail(z, A[i], design_tol, leveeFB)
     barrier_close <- dz*sum(fSLR[i:length(z)]*pf[i:length(z)]*damage[i:length(z)])
    
     aal_barrier[i] <- barrier_close + barrier_open
   }
   
   if (max(A) < target){
     print("Cannot find an allowance for barrier. The AAL target may be too large...")
     print(paste("AAL target=",target," max AAL levee=",max(aal_barrier),sep=""))
     return( NA )
   }
   return( A[which.min(abs(aal_barrier-target))] + leveeFB )
}

deci2lab <- function(x){
  if(x<1){
    y = x*10^(sum(floor(x*10^seq(0:16))!=x*10^seq(0:16))+1)
    return(paste("0p",y,sep=""))
  }else{
    return("1p0")  
  }
}

getCDF <- function(x,name){
  
  # smooth CDF
  dens = density(x, adjust=0.1, from=min(x) - 0.05*diff(range(x)), to=max(x) + 0.05*diff(range(x)))
  return( data.frame(CDF=cumsum(dens$y)/sum(dens$y),SLR=dens$x,name=name) )

  # empirical (step) CDF
  #f <- ecdf(x)
  #return( data.frame(CDF=f(x),SLR=x,AISmax=AISmax,name=name) )
}


Correct_Colnames <- function(df) {
  
  delete.columns <- grep("(^X$)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
  
  if (length(delete.columns) > 0) {
    
    row.names(df) <- as.character(df[, grep("^X$", colnames(df))])
    #other data types might apply than character or 
    #introduction of a new separate column might be suitable
    
    df <- df[,-delete.columns]
    
    colnames(df) <- gsub("^X", "",  colnames(df))
    #X might be replaced by different characters, instead of being deleted
  }
  
  return(df)
}

# Inundation depth-damage function

# Some depth-damage functions from fits to USACE

# North Atlanic Coast Comprehensive Study: Resilient Adaptation to Increasing Risk

# Physical depth-damage function summary report
# January 2015

# Commercial Building
combld <- function(h){
  if (h>0 & h<=3)
    return(min(0.0792+0.323*h+0.0526*h^2,1))
  else if (h>3)
    return(min(0.0792+0.323*3+0.0526*3^2,1))
  else
    return(0)
  end
}


# Two-Story Residential with Basement
res1bld <- function(h){
  if (h>0 & h<=3)
    return(min(0.18+0.178*h+0.0233*h^2-0.00778*h^3,1))
  else if (h>3)
    return(min(0.18+0.178*3+0.0233*3^2-0.00778*3^3,1))
  else
    return(0)
  end
}

# Urban Residential High Rise
res2bld <- function(h){
  if (h>0 & h<=3)
    return(min(0.142+0.0541*h-0.00368*h^2-0.00133*h^3,1))
  else if (h>3)
    return(min(0.142+0.0541*3-0.00368*3^2-0.00133*3^3,1))
  else
    return(0)
  end
}

depth_dmg <- function(z, bp, a1, b1, b2 ){

  # Inputs:
  # z is the extreme water level height (m above MHHW)

  # df is a data frame with:
  # 1. 1-D profile of accumulated property value
  # 2. Mid-points of the 1-D accumulated property value

  # Outputs:
  # value of damage for 1-D property profile given the extreme water level height, z

  dz <- diff(z)[1]
  zmax <- 10
  zz <- seq(0,zmax,dz)
  prop_margin1 <- sapply(zz,prop_margin,dz=dz,bp=bp,a1=a1,b1=b1,b2=b2) # marginal value of property
  
  nz <- length(z)
  ne <- length(z)
  
  dd <- sapply(z, function(x) res2bld(x)*0.95 + res1bld(x)*0.05)
  
  dmg_frac <- matrix(0, ne, nz)
  dmg_frac[lower.tri(dmg_frac, diag = TRUE)] <- dd[sequence(length(dd):1)]
  
  dmg_tot <- apply(sweep(dmg_frac,MARGIN=2,prop_margin1,`*`),1,sum)
  
  return(dmg_tot)
  }

# OLD Model damage vs "z" OLD
# depth_dmg <- function(z, df){
# 
#   # Inputs:
#   # z is the extreme water level height (m above MHHW)
#   
#   # df is a data frame with:
#   # 1. 1-D profile of accumulated property value
#   # 2. Mid-points of the 1-D accumulated property value
#   
#   # Outputs:
#   # value of damage for 1-D property profile given the extreme water level height, z
#   
#   de <- diff(df$e_list)[1]
#   
#   if (z > de) {
# 
#     dam1 <- NA
#     prop_margin <- diff(df$prop_acc) # marginal value of property
# 
#     i <- 1
#     for( e in seq(0,z,de) ){ # Get damage for this elevation (relative to MHHW)
#     
#       prop_e <- prop_margin[i] # property value of this elevation (mid-point)
#       h <- z - e # depth of the flood relative to elevation
#       print(h)
#       this_dam <- prop_e*(res2bld(h)*0.2 + res1bld(h)*0.5 + combld(h)*0.3)
# 
#       dam1[i] <- this_dam
#       #dam1 <- dam1 + this_dam
#       i <- i + 1
#     }
#     return(dam1)
#   }else{ # if z > 0
#     return(0)
#   } 
# }



depth_dmg_elev <- function(z, A, bp, a1, b1, a2, b2){
  
  # h_elev is the land elevation below which structures are elevated
  # A is the adaptation height (elevation height of the buildings)
  
  dz <- diff(z)[1]
  
 # # damage to elevated structures
 #  h <- max(0, z - A)
 #  dam1 <- prop_acc(A,bp,a1,b1,a2,b2)*(res2bld(h)*0.95 + res1bld(h)*0.05)
  
  dd <- sapply(pmax(0, z-A), function(x) res2bld(x)*0.95 + res1bld(x)*0.05)
  dam0 <- dd*prop_acc(A,bp,a1,b1,a2,b2)
  
  # damage to non-elevated structures
  dd <- sapply(z, function(x) res2bld(x)*0.95 + res1bld(x)*0.05)
  
  nz <- length(z)
  ne <- length(z)
  dmg_frac <- matrix(0, ne, nz)
  dmg_frac[lower.tri(dmg_frac, diag = TRUE)] <- dd[sequence(length(dd):1)]
  
#  zz <- c(0,rep(0,A/de),seq(A+de,zmax,dz))

#  dam1 <- apply(sweep(dmg_frac,MARGIN=2,prop_margin1,`*`),1,sum)
  
  prop_margin1 <- sapply(z,prop_margin,dz=dz,bp=bp,a1=a1,b1=b1,b2=b2) # marginal value of property
  tmp <- sweep(dmg_frac,MARGIN=2,prop_margin1,`*`)
  
  indx <- which.min(abs(z-A))
  dam1 <- apply(tmp[,indx:length(z)],1,sum)
 
  
  # sum damage to elevated and non-elevated
  return(dam1+dam0)

  # 
  # if (z > A){ # damage to non-elevated structures (above h)
  #   for(e in seq(A+de,z,de)){ 
  #       prop_e <- prop_margin(e,de,bp,a1,b1,b2) # value of property in this elevation bin
  #       h <- max(0, z - e)
  #       this_dam <- prop_e*(res2bld(h)*0.95 + res1bld(h)*0.05)
  #       dam1 <- dam1 + this_dam
  #   }
  # }
  # 
  # 
  # return(dam1)
}

prop_acc <- function( z, bp, a1, b1, a2, b2 ){
  # get modeled accumulated value of property for NYC (with break-point)
  # last updated: Wed Apr 24 21:26:24 EDT 2019
  
  if (z < bp){
    return( max(0,a1*max(0,z)+b1*max(0,z)^2) )
  }else{
    return( max(0, a2 + max(0,z)*b2) )
  }
}


prop_margin <- function( z, dz, bp, a1, b1, b2 ){
  # get modeled marginal value of property for NYC (with break-point)
  # last updated: Wed Apr 24 21:26:24 EDT 2019
  
  if (z <= bp){
    return( max(0,a1+2*b1*max(0,z))*dz )
  }else{
    return( b2*dz )
  }
}


D <- function( e, bp, a1, b1, a2, b2 ){
  # aggregate damage function for NYC (with break-point)
  # last updated: Wed Apr 24 21:26:24 EDT 2019
  
  if (z < bp){
    return( a1*max(0,z)^b1 )
    #return( max(0,a1*max(0,z)+b1*max(0,z)^2) )
  }else{
    return( max(0, a2 + max(0,z)*b2) )
  }
}