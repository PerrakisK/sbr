library(sbr)
library(mvtnorm)
library(pls)
library(normalp)

load('.../simulation_paper.RData') # load R enviroment which contains the covariance matrices of the simulated clinical and rna data
COV.RNA <- as.matrix(COV.RNA)
p.c <- 26
p.r <- 2000
p <- p.c + p.r                    # total number of predictors from the two sources (c: clinical, r: rna)
range.c <- 1:26
range.r <- 27:2026
nz.c <- p.c * 0.5
nz.r <- p.r * 0.05
set.seed(1); s.c <- sample(range.c, nz.c)
set.seed(1); s.r <- sample(range.r, nz.r)
beta.true <- rep(0, p)
beta.true[c(s.c, s.r)] <- 1
n <- 100                         # training sample size

# creating storage objects

repet <- 20                     # number of repetitions
Y <- matrix(NA, n, repet)
beta.enet <- matrix(NA, p, repet)
beta.lasso <- matrix(NA, p, repet)
beta.ssbr <- matrix(NA, p, repet)
beta.cssbr <- matrix(NA, p, repet)
beta.gssbr <- matrix(NA, p, repet)
grid.enet <- seq(0.1, 0.9, by = 0.1)
grid.length <- length(grid.enet)

# start of for loop - main calculations

for(i in 1:repet){
  
  
  set.seed(i); X.c <- rmvnorm(n, rep(0, p.c), COV.CL)
  set.seed(i); X.r <- rmvnorm(n, rep(0, p.r), COV.RNA)
  betas <- rep(0, p)
  sigma.c <- 0.5
  sigma.r <- 0.5
  betas[s.c] <- rnormp(nz.c, 0, sigma.c, p = 1.5)
  betas[s.r] <- rnormp(nz.r, 0, sigma.r, p = 1.5)
  set.seed(i); Y[, i] <-
    X.c[, s.c] %*% betas[s.c] + X.r[, s.r - 26] %*% betas[s.r] + rnorm(n, 0, 1)
  Y[, i] <- scale(Y[, i])
  X.c <- stdize(X.c)
  X.r <- stdize(X.r)
  X <- cbind(X.c, X.r)

  fold.id <- sample(rep(seq(5), length = n))
  model.enet <- list()
  for (l in 1:grid.length) {
    model.enet[[l]] <-                                 # elastic-net
      cv.glmnet(
        X,
        Y[, i],
        family = 'gaussian',
        standardize = FALSE,
        intercept = FALSE,
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
  beta.enet[, i] <-
    coef(model.enet[[ind.min.alpha]], s = 'lambda.min')[-1]
  
  model.lasso <-                                         # lasso
    cv.glmnet(
      X,
      Y[, i],
      nfolds = 5,
      family = 'gaussian',
      standardize = FALSE,
      intercept = FALSE
    )
  beta.lasso[, i] <- coef(model.lasso, s = 'lambda.min')[, 1][-1]
  
  Gram1 <- gram(X.c)
  Gram2 <- gram(X.r)
  X <- list(X.c, X.r)
  t_X.c <- t(X.c)
  t_X.r <- t(X.r)
  tX <- list(t_X.c, t_X.r)
  G <- list(Gram1, Gram2)
  
  model.ssbr <-                                           # relaxed SSBR
    sbr(
      Y[, i],
      X,
      tX,
      G,
      estimator = 'MAP',
      sparsify = TRUE,
      relaxed = TRUE,
      sparse.control = 1
    )
  beta.ssbr[, i] <- model.ssbr$coefficients[, 2]       
  
  model.cssbr <-                                          # relaxed cSSBR
    sbr(
      Y[, i],
      X,
      tX,
      G,
      estimator = 'MAP',
      sparsify = TRUE,
      relaxed = TRUE,
      sparse.control = log(n)
    )
  beta.cssbr[, i] <- model.cssbr$coefficients[, 2]
  
  model.gssbr <-                                          # general SSBR
    sbr(
      Y[, i],
      X,
      tX,
      G,
      estimator = 'MAP',
      sparsify = TRUE,
      relaxed = FALSE,
      sparse.control = 1
    )
  beta.gssbr[, i] <- model.gssbr$coefficients[, 2]
  
  print(i)
}

# exporting the betas

dput(beta.true, paste(path, 'true_', n, '_', '.txt', sep = ""))    # this are the true inclusion positions (binary vector)
dput(beta.enet, paste(path, 'enet_', n, '_', '.txt', sep = ""))
dput(beta.lasso, paste(path, 'lasso_', n, '_', '.txt', sep = ""))
dput(beta.ssbr, paste(path, 'ssbr_', n, '_', '.txt', sep = ""))
dput(beta.cssbr, paste(path, 'cssbr_', n, '_', '.txt', sep = ""))
dput(beta.gssbr, paste(path, 'gssbr_', n, '_', '.txt', sep = ""))



