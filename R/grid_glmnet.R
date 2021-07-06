# glmnet loop for fitting models
grid_glmnet = function( x, y, grid.size = p, foldassign = NULL, var_order = NULL, lambda = NULL, nlambda = 100, lambda.min.ratio = 0.01, K = 5, binary = F, maxit = 1e5, ...) {
  n = length(y)
  p = dim(x)[ 2 ]
  if ( is.null(foldassign) ) foldassign = ceiling(K * sample(1:n) / n)
  if ( is.null(var_order) ) var_order = order(apply(x, 2, var) )
  grid_seq = ( exp(log(p) / ( grid.size - 1 )) )^(seq.int(from = 0, to = grid.size - 1L, by = 1L))
  grid_seq = pmax(floor(grid_seq + 1e-8), rank(grid_seq, ties.method = "first")) # added rounding factor to ensure final term 
  # includes all variables
  grid_seq = p - grid_seq
  grid_seq = rev(grid_seq)
  # grid_seq = rev(rev(grid_seq)[ 1:grid.size.truncate ]) # Now that we have a cross-validation wrapper
  # grid.size = grid.size.truncate
  if ( grid.size == 1 ) {
    grid_seq[ 1 ] = 0
  }
  minerr = Inf
  for (m in 1:grid.size) {
    print(m)
    print(grid_seq[m])
    print(paste0("Excluding: ", length(head(var_order[ 1:grid_seq[ m ] ], grid_seq[ m ])), " variables."))
    # print(head(var_order[ 1:grid_seq[ m ] ], grid_seq[ m ]))
    newmodel = cv.glmnet(x, y, foldid = foldassign, thresh = 1e-10, family = ifelse(binary, "binomial", "gaussian"), 
                      lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, nfolds = K, parallel = T, exclude = head(var_order[ 1:grid_seq[ m ] ], grid_seq[ m ]), maxit = maxit)
    if ( min(newmodel$cvm) < minerr ) {
      minerr = min(newmodel$cvm)
      model = newmodel
      bestgrid = m
    }
  }
  model$gridpoint = bestgrid
  model$gridsize = p - grid_seq[ bestgrid ]
  return(model)
}

predict.grid_glmnet = function( object, newx, ... ) {
  return(predict(object, newx, s = "lambda.min", ...))
}