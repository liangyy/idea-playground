set.seed(2018)
source('em_algorithm.R')
source('simulate_data.R')
plot(abs(beta_eqtl), abs(beta_gwas))
sigma.e.em <- c()
sigma.1.em <- c()
sigma.2.em <- c()
pi.em <- c()
lld.em <- c()
for(i in 1 : 10){
  print(i)
  o <- em_algorithm(beta_gwas, beta_eqtl)
  sigma.e.em <- c(sigma.e.em, o$sigma)
  idx2 <- which.max(o$sigma_k)
  idx1 <- which.min(o$sigma_k)
  sigma.1.em <- c(sigma.1.em, o$sigma_k[idx1])
  sigma.2.em <- c(sigma.2.em, o$sigma_k[idx2])
  pi.em <- c(pi.em, o$pi[1])
  lld.em <- c(lld.em, max(o$lld))
}

df.em <- data.frame(sigma = sigma.e.em, sigma.1 = sigma.1.em, sigma.2 = sigma.2.em, pi = pi.em, lld = lld.em)
df.em$idx <- 1 : 10
lld <- get_loglik(beta_gwas, beta_eqtl, c(sigma.1, sigma.2), sigma.e, c(pi, 1 - pi))
df.true <- data.frame(sigma = sigma.e, sigma.1 = sigma.1, sigma.2 = sigma.2, pi = pi, lld = lld)

ggplot() + geom_jitter(data = melt(df.em, id.vars = 'idx'), aes(x = idx, y = value), width = 1, height = 0) + facet_wrap(~variable, scales = 'free') + geom_hline(data = melt(df.true), aes(yintercept = value), color = 'red')
