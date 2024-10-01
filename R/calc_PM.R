# file with functions to calculate performance metrics




#' compute Annual average variability in catch (AAV) see 
#'
#' @param ct vector of catches, bust be longer than 1
#' 
#' 
#' @returns Numeric AAV 
#' 
#' @details references: 
#' Cox, S.P. and Kronlund, A.R. 2008. Practical stakeholder-driven harvest policies 
#' for groundfish in British Columbia, Canada. Fish. Res. 94(3): 224-237
#' Punt, A.E., Smith, A.D.M., 1999. Harvest strategy evaluation for the eastern stock 
#' of gemfish (Rexea solnadri). ICES J. Mar. Sci. 56, 860–875
#' 
calc_AAV <- funcion(ct){

  if(length(ct<2)){stop("vector of catches 'ct' needs to be longer that 1")}

  AAV <- sum(abs(ct[1:(length(ct)-1)]-ct[2:(length(ct))]))/sum(ct)

  return(AAV)
}



