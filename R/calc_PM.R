# file with functions to calculate performance metrics


calc_AAV <- funcion(ct){

  if(length(ct<2)){stop("vector of catches 'ct' needs to be longer that 1")}

  AAV <- sum(abs(ct[1:(length(ct)-1)]-ct[2:(length(ct))]))/sum(ct)

  return(AAV)
}



