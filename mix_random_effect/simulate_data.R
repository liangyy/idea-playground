# simulate data

# beta_eqtl ~ N(0, sigma.beta_eqtl)
# beta_gene ~ pi N(0, sigma.1) + (1 - pi) N(0, sigma.2)
# e ~ N(0, sigma.e)

# n <- 10000
sigma.beta_eqtl <- 1
sigma.1 <- 0.001
sigma.2 <- 1
sigma.e <- 0.5
pi <- 0.7

beta_eqtl <- rnorm(n, mean = 0, sd = sigma.beta_eqtl)
z.gene <- runif(n) < pi
beta_gene <- rep(0, n)
beta_gene.1 <- rnorm(sum(z.gene), mean = 0, sd = sqrt(sigma.1))
beta_gene.2 <- rnorm(n - sum(z.gene), mean = 0, sd = sqrt(sigma.2))
beta_gene[z.gene] <- beta_gene.1
beta_gene[!z.gene] <- beta_gene.2
error <- rnorm(n, mean = 0, sd = sqrt(sigma.e))
beta_gwas <- beta_gene * beta_eqtl + error

