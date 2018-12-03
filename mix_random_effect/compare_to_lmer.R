# solving
# beta_gwas = r1 * beta_eqtl * I.1 + r2 * beta_eqtl * (1 - I.1) + error
# r1 ~ N(0, sigma.1)
# r2 ~ N(0, sigma.2)
# as sub-routine to solve 
# M-step optimization w.r.t sigma, sigma.1, sigma.2

# Proof of concept

set.seed(2018)

n <- 10000
source('simulate_data.R')
source('em_algorithm.R')

sigma_k <- c(sigma.1, sigma.2)
sigma <- sigma.e

beta_eqtl.square <- beta_eqtl ^ 2
logL <- sapply(sigma_k, function(sk) {dnorm(x = beta_gwas, mean = 0, sd = sqrt(beta_eqtl.square * sk + sigma), log = T)})
logw <- sweep(logL, 2, log(pis), FUN = '+')
logw.rowsum <- apply(logw, 1, logsum)
logw <- sweep(logw, 1, logw.rowsum, FUN = '-') 
logw <- matrix(log(0.5), ncol = 2, nrow = 10000)
out <- grad_solver(beta_gwas, beta_eqtl.square, logw, c(log(sigma_k), log(sigma)))
out2 <- exp(quad_solver(beta_gwas, beta_eqtl.square, logw, out$par))

library(lme4)
beta_eqtl.dummy <- cbind(c(beta_eqtl, rep(0, length(beta_eqtl))), c(rep(0, length(beta_eqtl)), beta_eqtl))
beta_gwas.dummy <- c(beta_gwas, beta_gwas)
df.lmer <- data.frame(y = beta_gwas.dummy, x1 = beta_eqtl.dummy[, 1], x2 = beta_eqtl.dummy[, 2], grp1 = 1 : length(beta_gwas.dummy), grp2 = 1 : length(beta_gwas.dummy), grp3 = 1 : length(beta_gwas.dummy))
out3 <- lmer(y ~ -1 + (0 + x1|grp1) + (0 + x2|grp2), df.lmer, REML = F, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE="ignore"))
out5 <- lmer(y ~ -1 + (0 + x1|grp1) + (0 + x2|grp2), df.lmer, REML = F, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE="ignore"), weights = c(exp(logw)[, 1], exp(logw)[,2]))
out4 <- as.data.frame(VarCorr(out3))
out6 <- as.data.frame(VarCorr(out5))
fn(log(out2), beta_gwas, beta_eqtl.square, logw)
fn(log(out4[,4]), beta_gwas, beta_eqtl.square, logw)
df.re <- data.frame(my_solver = out2, lmer = out4[, 4], lmer_weight_0.5 = out6[, 4])
row.names(df.re) <- c('sigma.1', 'sigma.2', 'sigma.e')
df.re

