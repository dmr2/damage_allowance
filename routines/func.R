ds <- function(s){
  
  return(0.894*s^2.23)
  
}

fs_ds <- function(s,dlt,mu,shp,scl,lmbd) {

  # Return GPD probability density function
  fs <- function(s){
    # This currently cannot handle surges off the support of GPD, 
    #  (e.g., s-mu-SLR < 0.51) which is a common issue with high values of SLR. 
    # 
    # A solution is needed...
    
    # y <- ifelse(s >= mu, gpd(s,mu,shp,scl),gumble(s,mu,scl))
    
    #if(any(s >= mu)){ # use the GPD PDF
    #if(s >= 0){ # use the GPD PDF
    exponent <- function(a, pow) (abs(a)^pow)*sign(a)
    return(lmbd/scl*exponent(1+(shp*pmax(s,0)/scl),(-1/shp)-1))
    
    # }else{ # use the Gumbel PDF
   # return(-lmbd/scl*exp(-(s-mu)/scl))
    
    # }
  }
  
  # Return damage for this amount of SLR
  ds <- function(s){ # meters only
    return(0.894*s^2.23) # return billions of USD$
  }
  
  return( ds(s)*fs(s-mu-dlt) )
} # end function FS_DS



getRSLsamps <- function(fil){
  
  # Open relative sea-level rise Monte Carlo samples (K14 format)
  
  # Output:
  #   - samples[ MC RSL samples (in meters), Years ]
  #   - Years
  
  
  print(" ")
  print( paste( "Opening SL MC samples:", fil) )
  print(" ")
  
  x = read.table( paste(fil, sep=""), skip=1, sep="\t", header=FALSE)
  
  years = x[,1] # Years corresponding to SLR
  
  iii <- as.data.frame( t(x) ) # Transpose to [SAMPLES,YEARS]
  iii <- iii[2:10001,]/1000 # millimeters to meters
  
  # remove samples that may be physically implausible (i.e., the 99.9th percentile)
  q99 <- apply( iii, 2, quantile, probs=c(0.999), na.rm=T)
  samples <- matrix( NaN, nrow=10000, ncol=length(years) )
  
  for(t in 1:length( years )){
    samples[,t] <- pmin(iii[,t],q99[t])
  }
  
  return ( list(samples=samples, years=years) )
  
}