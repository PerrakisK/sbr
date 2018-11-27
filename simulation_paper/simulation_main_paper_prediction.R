library(devtools)
devtools::install_github('kperrakis/sbr')
library(sbr)
library(mvtnorm)
library(pls)
library(normalp)
library(MCMCpack)

load('.../simulation_paper.RData') # load R enviroment which contains the covariance matrices of the simulated clinical and rna data
COV.RNA <- as.matrix(COV.RNA)
sigma.hat <- 0.1                 # error variance for the beta coefficients
p.c <- 26
p.r <- 2000
p.s <- 100000
p <- p.c + p.r + p.s            # total number of predictors from the three sources (c: clinical, r: rna, s: snp)
range.c <- 1:26
range.r <- 27:2026
range.s <- 2027:102026
nz.c <- p.c * 0.5
nz.r <- p.r * 0.05
nz1.s <- p.s * 0.01
nz2.s <- p.s * 0.1
nz3.s <- p.s * 0.5
set.seed(1); s.c <- sample(range.c, nz.c)
set.seed(1); s.r <- sample(range.r, nz.r)
set.seed(1); s1.s <- sample(range.s, nz1.s)   # sparsity levels scenario 1
set.seed(1); s2.s <- sample(range.s, nz2.s)   # sparsity levels scenario 2
set.seed(1); s3.s <- sample(range.s, nz3.s)   # sparsity levels scenario 3
block.size <- 100                             # block size in SNP covariance matrix
n.blocks <- 1000                              # number of blocks in SNP covariance matrix

# Two choices:
# 1) block.size=100 and n.blocks=1000 (low-correlation scenario)
# 2) block.size=1000 and n.blocks=100 (high correlation scenario)

if (n.blocks == 1000)
  correlation <- 'low'
if (n.blocks == 100)
  correlation <- 'high'

b.int <- matrix(NA, n.blocks, 2)
b.int[1, 1] <- range.s[1]
b.int[1, 2] <- range.s[1] + block.size - 1
for (j in 2:n.blocks) {
  b.int[j, 1] <- b.int[j - 1, 2] + 1
  b.int[j, 2] <- b.int[j, 1] + block.size - 1
}
COV.SNP <- list()
for (j in 1:n.blocks) {
  set.seed(j)
  COV.SNP[[j]] <-
    riwish(block.size, diag(block.size)) # block diagonal covariance matrix for the simulated SNP data
}
n <- 100                         # training sample size
n.test <- 5000                   # test sample size

# creating storage objects

repet <- 50
Y <- list()
Y.test <- list()
na.matrix <- matrix(NA, repet, 3)
colnames(na.matrix) <- c('sparse', 'medium', 'dense')
sparsity.ridge <- na.matrix
sparsity.lasso <- na.matrix
sparsity.enet <- na.matrix
sparsity.ml.ssbr <- na.matrix
sparsity.cv.ssbr <- na.matrix
sparsity.map.ssbr <- na.matrix
sparsity.ml.cssbr <- na.matrix
sparsity.cv.cssbr <- na.matrix
sparsity.map.cssbr <- na.matrix
cor.ridge <- na.matrix
cor.lasso <- na.matrix
cor.enet <- na.matrix
cor.ml <- na.matrix
cor.cv <- na.matrix
cor.map <- na.matrix
cor.ml.ssbr <- na.matrix
cor.cv.ssbr <- na.matrix
cor.map.ssbr <- na.matrix
cor.ml.cssbr <- na.matrix
cor.cv.cssbr <- na.matrix
cor.map.cssbr <- na.matrix
SD1 <- c()
SD2 <- c()
SD3 <- c()
for (j in 1:3) {
  Y[[j]] <- matrix(NA, n, repet)
  Y.test[[j]] <- matrix(NA, n.test, repet)
}
lambda.ml <- array(NA, c(repet, 3, 3))
lambda.cv <- array(NA, c(repet, 3, 3))
lambda.map <- array(NA, c(repet, 3, 3))
dimnames(lambda.ml)[[2]] <- c('CL', 'RNA', 'SNP')
dimnames(lambda.ml)[[3]] <- c('sparse', 'medium', 'dense')
dimnames(lambda.cv)[[2]] <- c('CL', 'RNA', 'SNP')
dimnames(lambda.cv)[[3]] <- c('sparse', 'medium', 'dense')
dimnames(lambda.map)[[2]] <- c('CL', 'RNA', 'SNP')
dimnames(lambda.map)[[3]] <- c('sparse', 'medium', 'dense')


# start of for loop - main calculations

for (i in 1:repet) {
  set.seed(i)
  X.c <- rmvnorm(n, rep(0, p.c), COV.CL)
  set.seed(i + 100)
  X.c.test <- rmvnorm(n.test, rep(0, p.c), COV.CL)
  set.seed(i)
  X.r <- rmvnorm(n, rep(0, p.r), COV.RNA)
  set.seed(i + 100)
  X.r.test <- rmvnorm(n.test, rep(0, p.r), COV.RNA)
  X.s <- list()
  X.s.test <- list()
  for (j in 1:n.blocks) {
    set.seed(j)
    X.s[[j]] <- rmvnorm(n, rep(0, block.size), COV.SNP[[j]])
    set.seed(j + 100)
    X.s.test[[j]] <-
      rmvnorm(n.test, rep(0, block.size), COV.SNP[[j]])
    X.s[[j]] <- abs(X.s[[j]])
    X.s[[j]][X.s[[j]] < 1.5] <- 0
    X.s[[j]][(X.s[[j]] >= 1.5) & (X.s[[j]] < 2.5)] <- 1
    X.s[[j]][X.s[[j]] >= 2.5] <- 2
    X.s.test[[j]] <- abs(X.s.test[[j]])
    X.s.test[[j]][X.s.test[[j]] < 1.5] <- 0
    X.s.test[[j]][(X.s.test[[j]] >= 1.5) &
                    (X.s.test[[j]] < 2.5)] <- 1
    X.s.test[[j]][X.s.test[[j]] >= 2.5] <- 2
  }
  X.s <- matrix(unlist(X.s), ncol = p.s, byrow = FALSE)
  X.s.test <- matrix(unlist(X.s.test), ncol = p.s, byrow = FALSE)
  betas <- matrix(0, p, 3)
  set.seed(i)
  betas[s.c, ] <- rnormp(nz.c, 0, sigma.hat, p = 1.5)
  betas[s.r, ] <- rnormp(nz.r, 0, sigma.hat, p = 1.5)
  betas[s1.s, 1] <- rnormp(nz1.s, 0, 2 * sigma.hat / 3, p = 1.5)
  betas[s2.s, 2] <- rnormp(nz2.s, 0, 2 * sigma.hat / 3, p = 1.5)
  betas[s3.s, 3] <- rnormp(nz3.s, 0, 2 * sigma.hat / 3, p = 1.5)
  set.seed(i); Y[[1]][, i] <-
    X.c[, s.c] %*% betas[s.c, 1] + X.r[, s.r - 26] %*% betas[s.r, 1] + X.s[, s1.s - 2026] %*% betas[s1.s, 1] + rnorm(n, 0, 1)
  set.seed(i); Y[[2]][, i] <-
    X.c[, s.c] %*% betas[s.c, 2] + X.r[, s.r - 26] %*% betas[s.r, 2] + X.s[, s2.s - 2026] %*% betas[s2.s, 2] + rnorm(n, 0, 1)
  set.seed(i); Y[[3]][, i] <-
    X.c[, s.c] %*% betas[s.c, 3] + X.r[, s.r - 26] %*% betas[s.r, 3] + X.s[, s3.s - 2026] %*% betas[s3.s, 3] + rnorm(n, 0, 1)
  set.seed(i + 100); Y.test[[1]][, i] <-
    X.c.test[, s.c] %*% betas[s.c, 1] + X.r.test[, s.r - 26] %*% betas[s.r, 1] + X.s.test[, s1.s - 2026] %*% betas[s1.s, 1] + rnorm(n, 0, 1)
  set.seed(i + 100); Y.test[[2]][, i] <-
    X.c.test[, s.c] %*% betas[s.c, 2] + X.r.test[, s.r - 26] %*% betas[s.r, 2] + X.s.test[, s2.s - 2026] %*% betas[s2.s, 2] + rnorm(n, 0, 1)
  set.seed(i + 100); Y.test[[3]][, i] <-
    X.c.test[, s.c] %*% betas[s.c, 3] + X.r.test[, s.r - 26] %*% betas[s.r, 3] + X.s.test[, s3.s - 2026] %*% betas[s3.s, 3] + rnorm(n, 0, 1)
  SD1[i] <- sd(Y[[1]][, i])
  SD2[i] <- sd(Y[[2]][, i])
  SD3[i] <- sd(Y[[3]][, i])
  Y[[1]][, i] <- (Y[[1]][, i] - mean(Y[[1]][, i])) / SD1[i]
  Y[[2]][, i] <- (Y[[2]][, i] - mean(Y[[2]][, i])) / SD2[i]
  Y[[3]][, i] <- (Y[[3]][, i] - mean(Y[[3]][, i])) / SD3[i]
  Y.test[[1]][, i] <-
    (Y.test[[1]][, i] - mean(Y.test[[1]][, i])) / sd(Y.test[[1]][, i])
  Y.test[[2]][, i] <-
    (Y.test[[2]][, i] - mean(Y.test[[2]][, i])) / sd(Y.test[[2]][, i])
  Y.test[[3]][, i] <-
    (Y.test[[3]][, i] - mean(Y.test[[3]][, i])) / sd(Y.test[[3]][, i])
  X.c <- stdize(X.c)
  X.r <- stdize(X.r)
  X.c.test <- stdize(X.c.test)
  X.r.test <- stdize(X.r.test)
  X.test <- cbind(X.c.test, X.r.test, X.s.test)
  t_X.c <- t(X.c)
  t_X.r <- t(X.r)
  t_X.s <- t(X.s)
  Gram1 <- gram(X.c)
  Gram2 <- gram(X.r)
  Gram3 <- gram(X.s)
  X <- list(X.c, X.r, X.s)
  tX <- list(t_X.c, t_X.r, t_X.s)
  G <- list(Gram1, Gram2, Gram3)
  XX <- cbind(X.c, X.r, X.s)
  grid.enet <- seq(0.1, 0.9, by = 0.1)
  grid.length <- length(grid.enet)
  fold.id <- sample(rep(seq(5), length = n))
  model.enet <- list()
  
  for (j in 1:3) {
    
    n.cores <- 5                 # number of cores for 5-fold CV
    cl <- makeCluster(n.cores)   # make cluster for running elastic-net, ridge, lasso
    registerDoParallel(cl)
    
    for (l in 1:grid.length) {
      model.enet[[l]] <-
        cv.glmnet(
          XX,
          Y[[j]][, i],
          family = 'gaussian',
          standardize = FALSE,
          intercept = FALSE,
          parallel = TRUE,
          alpha = grid.enet[l],
          foldid = fold.id
        )
    }
    ind.cv <-
      lapply(1:grid.length, function(x)
        which.min(model.enet[[x]]$cvm))
    ind.min.alpha <-
      which.min(lapply(1:grid.length, function(x)
        model.enet[[x]]$cvm[ind.cv[[x]]]))
    beta.enet <-
      coef(model.enet[[ind.min.alpha]], s = 'lambda.min')[-1]
    sparsity.enet[i, j] <-
      length(which(beta.enet != 0)) / p
    pred.enet <- X.test %*% beta.enet
    cor.enet[i, j] <-
      cor(Y.test[[j]][, i], pred.enet, use = 'complete.obs')
    
    model.ridge <-
      cv.glmnet(
        XX,
        Y[[j]][, i],
        nfolds = 5,
        family = 'gaussian',
        standardize = FALSE,
        intercept = FALSE,
        parallel = TRUE,
        alpha = 0
      )
    beta.ridge <- coef(model.ridge, s = 'lambda.min')[, 1][-1]
    pred.ridge <- X.test %*% beta.ridge
    cor.ridge[i, j] <-
      cor(Y.test[[j]][, i], pred.ridge, use = 'complete.obs')
    
    model.lasso <-
      cv.glmnet(
        XX,
        Y[[j]][, i],
        nfolds = 5,
        family = 'gaussian',
        standardize = FALSE,
        parallel = TRUE,
        intercept = FALSE
      )
    beta.lasso <- coef(model.lasso, s = 'lambda.min')[, 1][-1]
    sparsity.lasso[i, j] <-
      length(which(beta.lasso != 0)) / p
    pred.lasso <- X.test %*% beta.lasso
    cor.lasso[i, j] <-
      cor(Y.test[[j]][, i], pred.lasso, use = 'complete.obs')
    
    stopCluster(cl)  # stopping cluster, continuing with SBR/SSBR methods
    
    model.ml.ssbr <-
      sbr(
        Y[[j]][, i],
        X,
        tX,
        G,
        estimator = 'ML',
        sparsify = TRUE,
        relaxed = TRUE,
        sparse.control = 1
      )
    model.ml.cssbr <-
      sbr(
        Y[[j]][, i],
        X,
        tX,
        G,
        estimator = 'ML',
        sparsify = TRUE,
        relaxed = TRUE,
        sparse.control = log(n)
      )
    model.cv.ssbr <-
      sbr(
        Y[[j]][, i],
        X,
        tX,
        G,
        estimator = 'CV',
        sparsify = TRUE,
        relaxed = TRUE,
        sparse.control = 1
      )
    model.cv.cssbr <-
      sbr(
        Y[[j]][, i],
        X,
        tX,
        G,
        estimator = 'CV',
        sparsify = TRUE,
        relaxed = TRUE,
        sparse.control = log(n)
      )
    model.map.ssbr <-
      sbr(
        Y[[j]][, i],
        X,
        tX,
        G,
        estimator = 'MAP',
        sparsify = TRUE,
        relaxed = TRUE,
        sparse.control = 1
      )
    model.map.cssbr <-
      sbr(
        Y[[j]][, i],
        X,
        tX,
        G,
        estimator = 'MAP',
        sparsify = TRUE,
        relaxed = TRUE,
        sparse.control = log(n)
      )
    
    beta.ml <- coef(model.ml.ssbr)[, 1]
    beta.ml.ssbr <- coef(model.ml.ssbr)[, 2]
    beta.ml.cssbr <- coef(model.ml.cssbr)[, 2]
    sparsity.ml.ssbr[i, j] <- length(which(beta.ml.ssbr != 0)) / p
    sparsity.ml.cssbr[i, j] <- length(which(beta.ml.cssbr != 0)) / p
    pred.ml <- X.test %*% beta.ml
    pred.ml.ssbr <- X.test %*% beta.ml.ssbr
    pred.ml.cssbr <- X.test %*% beta.ml.cssbr
    cor.ml[i, j] <- cor(Y.test[[j]][, i], pred.ml)
    cor.ml.ssbr[i, j] <- cor(Y.test[[j]][, i], pred.ml.ssbr)
    cor.ml.cssbr[i, j] <- cor(Y.test[[j]][, i], pred.ml.cssbr)
    
    beta.cv <- coef(model.cv.ssbr)[, 1]
    beta.cv.ssbr <- coef(model.cv.ssbr)[, 2]
    beta.cv.cssbr <- coef(model.cv.cssbr)[, 2]
    sparsity.cv.ssbr[i, j] <- length(which(beta.cv.ssbr != 0)) / p
    sparsity.cv.cssbr[i, j] <- length(which(beta.cv.cssbr != 0)) / p
    pred.cv <- X.test %*% beta.cv
    pred.cv.ssbr <- X.test %*% beta.cv.ssbr
    pred.cv.cssbr <- X.test %*% beta.cv.cssbr
    cor.cv[i, j] <- cor(Y.test[[j]][, i], pred.cv)
    cor.cv.ssbr[i, j] <- cor(Y.test[[j]][, i], pred.cv.ssbr)
    cor.cv.cssbr[i, j] <- cor(Y.test[[j]][, i], pred.cv.cssbr)
    
    beta.map <- coef(model.map.ssbr)[, 1]
    beta.map.ssbr <- coef(model.map.ssbr)[, 2]
    beta.map.cssbr <- coef(model.map.cssbr)[, 2]
    sparsity.map.ssbr[i, j] <- length(which(beta.map.ssbr != 0)) / p
    sparsity.map.cssbr[i, j] <- length(which(beta.map.cssbr != 0)) / p
    pred.map <- X.test %*% beta.map
    pred.map.ssbr <- X.test %*% beta.map.ssbr
    pred.map.cssbr <- X.test %*% beta.map.cssbr
    cor.map[i, j] <- cor(Y.test[[j]][, i], pred.map)
    cor.map.ssbr[i, j] <- cor(Y.test[[j]][, i], pred.map.ssbr)
    cor.map.cssbr[i, j] <- cor(Y.test[[j]][, i], pred.map.cssbr)
    
    lambda.ml[i, , j] <- model.ml.ssbr$lambda
    lambda.cv[i, , j] <- model.cv.ssbr$lambda
    lambda.map[i, , j] <- model.map.ssbr$lambda
  }
  print(i)
}

path <-
  'C:/Users/user-adm/Desktop/simulation_paper/'  # path to save output

# sparsity levels

dput(sparsity.enet,
     paste(path, 'sparsity_enet_', n, '_', correlation, '.txt', sep = ""))
dput(sparsity.lasso,
     paste(path, 'sparsity_lasso_', n, '_', correlation, '.txt', sep = ""))
dput(
  sparsity.ml.ssbr,
  paste(path, 'sparsity_ml_ssbr_', n, '_', correlation, '.txt', sep = "")
)
dput(
  sparsity.cv.ssbr,
  paste(path, 'sparsity_cv_ssbr_', n, '_', correlation, '.txt', sep = "")
)
dput(
  sparsity.map.ssbr,
  paste(path, 'sparsity_map_ssbr_', n, '_', correlation, '.txt', sep = "")
)
dput(
  sparsity.ml.cssbr,
  paste(path, 'sparsity_ml_cssbr_', n, '_', correlation, '.txt', sep = "")
)
dput(
  sparsity.cv.cssbr,
  paste(path, 'sparsity_cv_cssbr_', n, '_', correlation, '.txt', sep = "")
)
dput(
  sparsity.map.cssbr,
  paste(path, 'sparsity_map_cssbr_', n, '_', correlation, '.txt', sep = "")
)

# correlations

dput(cor.enet,
     paste(path, 'cor_enet_', n, '_', correlation, '.txt', sep = ""))
dput(cor.lasso,
     paste(path, 'cor_lasso_', n, '_', correlation, '.txt', sep = ""))
dput(cor.ridge,
     paste(path, 'cor_ridge_', n, '_', correlation, '.txt', sep = ""))
dput(cor.ml,
     paste(path, 'cor_ml_', n, '_', correlation, '.txt', sep = ""))
dput(cor.cv,
     paste(path, 'cor_cv_', n, '_', correlation, '.txt', sep = ""))
dput(cor.map,
     paste(path, 'cor_map_', n, '_', correlation, '.txt', sep = ""))
dput(cor.ml.ssbr,
     paste(path, 'cor_ml_ssbr_', n, '_', correlation, '.txt', sep = ""))
dput(cor.cv.ssbr,
     paste(path, 'cor_cv_ssbr_', n, '_', correlation, '.txt', sep = ""))
dput(cor.map.ssbr,
     paste(path, 'cor_map_ssbr_', n, '_', correlation, '.txt', sep = ""))
dput(cor.ml.cssbr,
     paste(path, 'cor_ml_cssbr_', n, '_', correlation, '.txt', sep = ""))
dput(cor.cv.cssbr,
     paste(path, 'cor_cv_cssbr_', n, '_', correlation, '.txt', sep = ""))
dput(cor.map.ssbr,
     paste(path, 'cor_map_cssbr_', n, '_', correlation, '.txt', sep = ""))

# SBR lambda estimates

dput(lambda.ml,
     paste(path, 'lambda_ml_', n, '_', correlation, '.txt', sep = ""))
dput(lambda.cv,
     paste(path, 'lambda_cv_', n, '_', correlation, '.txt', sep = ""))
dput(lambda.map,
     paste(path, 'lambda_map_', n, '_', correlation, '.txt', sep = ""))
