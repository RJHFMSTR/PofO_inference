# Internal Function
#
# Calculate IBD coefficients.
# @param identity.coefficients a matrix with 1 row and 9 columns.
# @return The IBD coefficients.
ibd_coefs_A <- function(identity.coefficients){
  z0 <- identity.coefficients[,4] + identity.coefficients[,6] + identity.coefficients[,8] + identity.coefficients[,11]
  z1 <- identity.coefficients[,5] + identity.coefficients[,7] + identity.coefficients[,10]
  z2 <- identity.coefficients[,3] + identity.coefficients[,9]
  return(cbind(z0, z1, z2))
}
