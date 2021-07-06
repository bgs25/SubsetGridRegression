#This code is lifted from BDCoCoLasso package

ADMM_proj = function (mat, epsilon = 1e-04, mu = 10, it.max = 1000, etol = 1e-04, 
          etol_distance = 1e-04) 
{
  p <- nrow(mat)
  R <- diag(mat)
  S <- matrix(0, p, p)
  L <- matrix(0, p, p)
  itr <- 0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr < it.max) {
    Rp <- R
    Sp <- S
    start <- Sys.time()
    W <- mat + S + mu * L
    W.eigdec <- eigen(W, symmetric = TRUE)
    W.V <- W.eigdec$vectors
    W.D <- W.eigdec$values
    R <- W.V %*% diag(pmax(W.D, epsilon)) %*% t(W.V)
    M <- R - mat - mu * L
    S[lower.tri(S, diag = TRUE)] <- M[lower.tri(M, diag = TRUE)] - 
      l1proj(v = M[lower.tri(M, diag = TRUE)], b = mu/2)
    for (i in 2:p) {
      for (j in 1:(i - 1)) {
        S[j, i] <- S[i, j]
      }
    }
    L <- L - (R - S - mat)/mu
    end <- Sys.time()
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R, max(abs(R - Rp)))
    eps_S <- c(eps_S, max(abs(S - Sp)))
    eps_primal <- c(eps_primal, max(abs(R - S - mat)))
    time <- c(time, end - start)
    distance <- c(distance, max(abs(R - mat)))
    if (((max(abs(R - Rp)) < etol) && (max(abs(S - Sp)) < 
                                       etol) && (max(abs(R - S - mat)) < etol)) || (abs(max(abs(Rp - 
                                                                                                mat)) - max(abs(R - mat))) < etol_distance)) {
      itr <- it.max
    }
    else {
      itr <- itr + 1
    }
    if (itr%%20 == 0) {
      mu <- mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, 
                        eps_S = eps_S, eps_primal = eps_primal, time = time, 
                        distance = distance)
  return(list(mat = R, df_ADMM = df_ADMM))
}

l1proj = function (v, b) 
{
  stopifnot(b > 0)
  u <- sort(abs(v), decreasing = TRUE)
  sv <- cumsum(u)
  rho <- max(which(u > (sv - b)/1:length(u)))
  theta <- max(0, (sv[rho] - b)/rho)
  w <- sign(v) * pmax(abs(v) - theta, 0)
  return(w)
}



