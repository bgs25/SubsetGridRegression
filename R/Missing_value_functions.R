matvecprod = function( x, y ) {
  return(sapply(1:dim(x)[ 1 ], function(j) mean(x[ j, ] * y, na.rm = T)))
}

# matvecprod3 = function( x, y ) {
#   return(sapply(1:dim(x)[ 1 ], function(j) mean3(x[ j, ] * y, na.rm = T)))
# }

# matmatprod = function( x1, x2 ) {
#   return(sapply(1:dim(x2)[ 2 ], function(j) matvecprod( x1, x2[ , j ])))
# }

# matmatprod2 = function( x1, x2 ) {
#   return(sapply(1:dim(x2)[ 2 ], function( j ) sapply(1:dim(x1)[ 1 ], function(k) mean(x1[ k, ] * x2[ , j ], na.rm = T))))
# }
        
# matmatprod3 = function( x1 ) {
#   x1miss = !is.na(x1)
#   x1[ is.na(x1) ] = 0
#   return(crossprod(x1) / (crossprod(x1miss) - 1))
# }

matmatprod = function( x1 ) {
  x1miss = !is.na(x1)
  x1[ is.na(x1) ] = 0
  return(base::crossprod(x1) / (crossprod(x1miss)))
}

# mean3 = function( x, na.rm = T ) {
#   return(sum(x, na.rm = T) / (sum(!is.na(x))-1))
# }

scale_mle = function( x ) {
  return(scale(x, scale = sqrt(apply(scale(x,scale=F)^2, 2, mean, na.rm = T))))
}