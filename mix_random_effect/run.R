set.seed(2018)
library(ggplot2)
library(reshape2)
library(gridExtra)

source('em_algorithm.R')

plist <- list()

for(n in 3 : 5) {
  n <- 10 ^ n
  source('simulate_data.R')
  plist[[length(plist) + 1]] <- ggplot(data.frame(beta_eqtl = beta_eqtl, beta_gwas = beta_gwas)) + geom_point(aes(x = abs(beta_eqtl), y = abs(beta_gwas))) + ggtitle(paste0('n = ', n))
  sigma.e.em <- c()
  sigma.1.em <- c()
  sigma.2.em <- c()
  pi.em <- c()
  lld.em <- c()
  for(i in 1 : 10){
    print(i)
    o <-  tryCatch({
      em_algorithm(beta_gwas, beta_eqtl, seed = i)
    }, error = function(e) {
      return(NULL)
    })
    if(!is.null(o)) {
      sigma.e.em <- c(sigma.e.em, o$sigma)
      idx2 <- which.max(o$sigma_k)
      idx1 <- which.min(o$sigma_k)
      sigma.1.em <- c(sigma.1.em, o$sigma_k[idx1])
      sigma.2.em <- c(sigma.2.em, o$sigma_k[idx2])
      pi.em <- c(pi.em, o$pi[idx1])
      lld.em <- c(lld.em, max(o$lld))
    }
  }
  
  df.em <- data.frame(sigma = sigma.e.em, sigma.1 = sigma.1.em, sigma.2 = sigma.2.em, pi = pi.em, lld = lld.em)
  df.em$idx <- 1 : length(sigma.e.em)
  lld <- get_loglik(beta_gwas, beta_eqtl, c(sigma.1, sigma.2), sigma.e, c(pi, 1 - pi))
  df.true <- data.frame(sigma = sigma.e, sigma.1 = sigma.1, sigma.2 = sigma.2, pi = pi, lld = lld)
  
  plist[[length(plist) + 1]] <- ggplot() + geom_point(data = melt(df.em, id.vars = 'idx'), aes(x = factor(idx), y = value)) + facet_wrap(~variable, scales = 'free', nrow = 1) + geom_hline(data = melt(df.true), aes(yintercept = value), color = 'red') + ggtitle(paste0('n = ', n))
}

do.call("grid.arrange", c(plist[1:2], ncol = 1))
do.call("grid.arrange", c(plist[3:4], ncol = 1))
do.call("grid.arrange", c(plist[5:6], ncol = 1))