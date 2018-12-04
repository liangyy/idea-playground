em_algorithm <- function(beta_gwas, beta_eqtl, K = 2, tol = 1e-5, init = NULL, max.niter = 100, seed = 1) {
  set.seed(seed)
  niter <- 0
  # initialization
  if(is.null(init)) {
    sigma_k <- runif(K)
    sigma <- runif(1)
    pi <- runif(K)
    pi <- pi / sum(pi)
  } else {
    sigma_k <- init$sigma_k
    sigma <- init$sigma
    pi <- init$pi
  }
  
  loglik <- -Inf
  lld <- c()
  beta_eqtl.square <- beta_eqtl ^ 2
  loglik.max <- loglik
  
  while(max.niter > niter) {
    # E-step
    
    logL <- sapply(sigma_k, function(sk) {dnorm(x = beta_gwas, mean = 0, sd = sqrt(beta_eqtl.square * sk + sigma), log = T)})
    logw <- sweep(logL, 2, log(pi), FUN = '+')
    logw.rowsum <- apply(logw, 1, logsum)
    logw <- sweep(logw, 1, logw.rowsum, FUN = '-') 
    loglik.new.element <- sweep(logL, 2, log(pi), FUN = '+')
    loglik.new.rowsum <- apply(loglik.new.element, 1, logsum)
    loglik.new <- sum(loglik.new.rowsum)
    lld <- c(lld, loglik.new)
    # message('loglik = ', loglik, ', loglik.new = ', loglik.new)
    # message('sigma.k = ', paste(sigma_k, collapse = ', '), ' sigma = ', sigma, ' pi = ', paste(pi, collapse = ', '))
    if(loglik.new > loglik.max) {
      loglik.max <- loglik.new
      sigma.out = sigma
      sigma_k.out = sigma_k
      pi.out = pi
    }
    if(loglik.new < loglik) {
      # message('Error: does not converge!')
      loglik <- loglik.new
      # return(list(sigma = sigma.o, sigma_k = sigma_k.o, pi = pi.o, lld = lld))
    } else if(abs(loglik - loglik.new) > tol) {
      loglik <- loglik.new
    } else {
      break
    }
    
    # M-step
    logw.colsum <- apply(logw, 2, logsum)
    # message('dim logw.colsum = ', paste0(dim(logw.colsum), collapse = ' '))
    log2.sum <- logsum(logw.colsum)
    # message('dim log2.sum = ', paste0(dim(log2.sum), collapse = ' '))
    lld.before <- get_loglik(beta_gwas, beta_eqtl, sigma_k, sigma, pi)
    # q1 <- sum(sweep(exp(logw), 2, log(pi), FUN = '*') )
    pi <- exp(logw.colsum - log2.sum)
    # q2 <- sum(sweep(exp(logw), 2, log(pi), FUN = '*') )
    # if(q2 - q1 < 0) {
      
      # message('del q = ', q2 - q1)
    # }
    # lld.a1 <- get_loglik(beta_gwas, beta_eqtl, sigma_k, sigma, pi)
    # fn.before <- fn(c(log(sigma_k), log(sigma)), beta_gwas, beta_eqtl.square, logw)
    
    out <- grad_solver(beta_gwas, beta_eqtl.square, logw, c(log(sigma_k), log(sigma)))
    # message('grad_fn (optim) = ', paste(grad_fn(out$par, beta_gwas, beta_eqtl.square, logw), collapse = ', '))
    out2 <- quad_solver(beta_gwas, beta_eqtl.square, logw, out$par)
    # message('grad_fn (newton) = ', paste(grad_fn(out2, beta_gwas, beta_eqtl.square, logw), collapse = ', '))
    # fn.after <- fn(out2, beta_gwas, beta_eqtl.square, logw)
    # message('grad_fn (newton) = ', paste(grad_fn(out2, beta_gwas, beta_eqtl.square, logw), collapse = ', '))
    # message('del fn = ', fn.after - fn.before)
    sigma_k <- exp(out2[1 : (length(out$par) - 1)])
    sigma <- exp(out2[length(out$par)])
    # lld.a2 <- get_loglik(beta_gwas, beta_eqtl, sigma_k, sigma, pi)
    # message('del lld1 = ', lld.a1 - lld.before)
    # message('del lld2 = ', lld.a2 - lld.a1)
    niter <- niter + 1
    # print(loglik)
  }
  return(list(sigma = sigma.out, sigma_k = sigma_k.out, pi = pi.out, lld = lld, posterior_z = get_posterior_z(beta_gwas, beta_eqtl, sigma_k.out, sigma.out, pi.out)))
}

logsum <- function(logx) {
  # return log(sum(x))
  max.log.x <- max(logx)
  x <- sum(exp(logx - max.log.x))
  return(log(x) + max.log.x)
}

get_tau <- function(beta_eqtl.square, sigma_k, sigma){
  o <- matrix(beta_eqtl.square, ncol = 1, nrow = length(beta_eqtl.square)) %*% matrix(sigma_k, ncol = length(sigma_k), nrow = 1) + sigma
  return(o)
}

grad_fn <- function(params, beta_gwas, beta_eqtl.square, logw) {
  x_k <- params[1 : (length(params) - 1)]
  x <- params[length(params)]
  tau <- get_tau(beta_eqtl.square, exp(x_k), exp(x))
  w <- exp(logw)
  betaw <- sweep(w, 1, beta_eqtl.square, FUN = '*')
  factor <- -1 / tau + sweep(1 / tau ^ 2, 1, beta_gwas ^ 2, FUN = '*')
  grad_k <- colSums(betaw * factor) * 0.5 * exp(x_k)
  grad_s <- sum(w * factor) * 0.5 * exp(x)
  return(c(grad_k, grad_s))
}

fn <- function(params, beta_gwas, beta_eqtl.square, logw) {
  x_k <- params[1 : (length(params) - 1)]
  x <- params[length(params)]
  tau <- get_tau(beta_eqtl.square, exp(x_k), exp(x))
  w <- exp(logw)
  o <- sum(0.5 * w * (log(1 / tau) - sweep(1 / tau, 1, beta_gwas ^ 2, FUN = '*')))
  return(o)
}

grad_solver <- function(beta_gwas, beta_eqtl.square, logw, params) {
  optim(params, fn = function(par) { fn(par, beta_gwas, beta_eqtl.square, logw)}, method = "L-BFGS-B", gr = function(par) { grad_fn(par, beta_gwas, beta_eqtl.square, logw)}, control = list(fnscale = -1), hessian = T) # 
}

quad_solver <- function(beta_gwas, beta_eqtl.square, logw, params) {
  local_root_finding(x.init = params, fn = function(par) { grad_fn(par, beta_gwas, beta_eqtl.square, logw)}, jn = function(par) {jacobian_grad(par, beta_gwas, beta_eqtl.square, logw)})
}

local_root_finding <- function(x.init, fn, jn, tol = 1e-10, max.niter = 10) {
  diff <- Inf
  x.n <- x.init
  niter <- 1
  l.fn <- sqrt(sum(fn(x.n)^2))
  while(max.niter > niter) {
    del.x = tryCatch({
      solve(jn(x.n), fn(x.n))
    }, error = function(e) {
      return(x.n)
    })
    x.new <- x.n - del.x
    l.fn.new <- sqrt(sum(fn(x.new)^2))
    if(l.fn.new > l.fn) return(x.n)
    # message(paste(del.x, collapse = ' '))
    # message(paste(fn(x.new), collapse = ' '))
    diff <- sqrt(sum((x.n - x.new)^2))
    if(diff == Inf) return(x.n)
    if(diff < tol) {
      return(x.new)
    } else {
      x.n <- x.new
    }
    niter <- niter + 1
    l.fn <- l.fn.new
  }
  return(x.new)
}

jacobian_grad <- function(params, beta_gwas, beta_eqtl.square, logw) {
  x_k <-params[1 : (length(params) - 1)]
  x <- params[length(params)]
  tau <- get_tau(beta_eqtl.square, exp(x_k), exp(x))
  w <- exp(logw)
  betaw <- sweep(w, 1, beta_eqtl.square, FUN = '*')
  beta2w <- sweep(w, 1, beta_eqtl.square ^ 2, FUN = '*')
  factor <- - 1 / tau + sweep(1 / tau ^ 2, 1, beta_gwas ^ 2, FUN = '*')
  factor2 <- 1 / tau ^ 2 - sweep(1 / tau ^ 3, 1, 2 * beta_gwas ^ 2, FUN = '*')
  grad_kk <- colSums(betaw * factor) * 0.5 * exp(x_k) + colSums(beta2w * factor2) * 0.5 * exp(2 * x_k)
  grad_ks <- colSums(betaw * factor2) * 0.5 * exp(x_k) * exp(x)
  grad_ss <- sum(w * factor) * 0.5 * exp(x) + sum(w * factor2) * 0.5 * exp(2 * x)
  j_kk <- diag(grad_kk)
  j <- cbind(rbind(j_kk, grad_ks), c(grad_ks, grad_ss))
  return(j)
}

get_loglik <- function(beta_gwas, beta_eqtl, sigma_k, sigma, pi) {
  beta_eqtl.square <- beta_eqtl ^ 2
  logL <- sapply(sigma_k, function(sk) {dnorm(x = beta_gwas, mean = 0, sd = sqrt(beta_eqtl.square * sk + sigma), log = T)})
  loglik.new.element <- sweep(logL, 2, log(pi), FUN = '+')
  loglik.new.rowsum <- apply(loglik.new.element, 1, logsum)
  loglik.new <- sum(loglik.new.rowsum)
  return(loglik.new)
}

get_posterior_z <- function(beta_gwas, beta_eqtl, sigma_k, sigma, pi) {
  beta_eqtl.square <- beta_eqtl ^ 2
  logL <- sapply(sigma_k, function(sk) {dnorm(x = beta_gwas, mean = 0, sd = sqrt(beta_eqtl.square * sk + sigma), log = T)})
  logw <- sweep(logL, 2, log(pi), FUN = '+')
  logw.rowsum <- apply(logw, 1, logsum)
  logw <- sweep(logw, 1, logw.rowsum, FUN = '-') 
  return(exp(logw)[, 1])
}
