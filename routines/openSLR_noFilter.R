getRSLsamps_noFilter <- function(fil){
  
  # Open relative sea-level rise Monte Carlo samples (K14 format)
  
  # Output:
  #   - samples MC RSL samples (in meters), Years ]
  #   - Years

  
  print(" ")
  print( paste( "Opening SL MC samples:", fil) )
  print(" ")
  
  
  x = read.table( fil, skip=1, sep="\t", header=FALSE)
  
  years = x[,1] # Years corresponding to SLR
  
  iii <- as.data.frame( t(x) ) # Transpose to [SAMPLES,YEARS]
  samples <- iii[2:10001,]/1000 # millimeters to meters
  
  return ( list(samples=samples, years=years) )
  
}