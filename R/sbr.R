#' Scalable Bayesian regression (SBR) and sparse SBR for normal linear regression with multiple predictor matrices.
#' 
#' Function for scalable Bayesian regression (SBR/SSBR) in normal linear models with multiple types (sources) of feature matrices (with K being the number of sources). When K = 1, SBR corresponds to standard ridge regression using one from the three available empirical Bayes estimators (see below) for the penalty parameter. For details see
#' Perrakis and Mukherjee (2018).
#' @author Konstanstinos Perrakis \email{konstantinos.perrakis@dzne.de}
#' @author Sach Mukherjee \email{sach.mukherjee@dzne.de}
#' @references Perrakis, K., Mukherjee, S and the Alzheimers Disease Neuroimaging Initiative. (2018) Scalable Bayesian regression in high dimensions with multiple data sources. \url{https://arxiv.org/pdf/1710.00596.pdf}
#' @seealso \code{\link{gram}}
#' @seealso \code{\link{gram.parallel}}
#' @seealso \code{\link{predict.sbr}}
#' @param y a standardized response vector. 
#' @param X a standardized feature matrix (if K = 1) or a list of standardized feature matrices (if K > 1).
#' @param trX (optional) the transpose matrix of X (if K = 1) or a list of transpose matrices (if K > 1).
#' @param G the inner-product Gram matrix (if K = 1) or a list containing the multiple Gram matrices (if K > 1).
#' @param estimator the estimator used for tuning the shrinkage levels. Available estimates are leave-one-out cross-validation ("CV"), maximum marginal likelihood ("ML") and the maximum-a-posteriori value ("MAP", default).
#' @param sparsify logical, if TRUE the SSBR solution is calculated, default option is FALSE.
#' @param relaxed logical, if TRUE (default) the relaxed SSBR solution is calculated, if FALSE the general SSBR solution is calculated.
#' @param sparse.control numerical value for controlling the effect of sample size (n) on the resulting SSBR solution. Default option is 1 (no control). A recommended option for sparser solutions is sparse.control = log(n).
#' @param p.threshold used for block-matrix computation of the main diagonal of the covariance matrix when sparsify = TRUE and relaxed = TRUE. It will be triggered for any source-matrix whose number of columns is larger than p.threshold.   
#' @param cov.blocks argument corresponding to block size (not the number of blocks) when the block-matrix computation is triggered (see above). Default option is 1000, i.e. blocks of dimensionality 1000 x 1000.
#' @param parallel logical, if parallel = TRUE the calculation of variance components required for the SSBR solution is performed in parallel. Default is FALSE.
#' @param cl the number of cores to use when parallel = TRUE. Must be provided by the user.
#' @param L.optim lower bound for the optimization procedure used to tune the shrinkage levels, default is 1e-04.
#' @param U.optim upper bound for the optimization procedure used to tune the shrinkage levels, default is 1e04.
#' @return An object with S3 class 'sbr' allowing for call to generic functions \code{\link{coef}} and \code{\link{predict}}  
#' @return An object of class 'sbr' is a list containing the following components:
#' \item{coefficients}{a 1-column matrix with the SBR beta estimates (when sparsify = FALSE) or a 2-column matrix with the SBR and SSBR beta estimates (when sparsify = TRUE). Note that the coefficients correspond to the standardized response variable and feature matrix.}
#' \item{sigma2}{the variance component (at the posterior mode).} 
#' \item{lambda}{the vector of penalty parameters.}
#' \item{lambdaEstimator}{the estimation method for lambda.}
#' \item{duration}{reported runtime}
#' @export 
#' @examples 
#' ################# Toy example with 3 data sources #################
#' library(mvtnorm)
#' library(MCMCpack)
#' ### GENERATION OF DATA ###
#' ## sample size and number of predictors
#' n<-50 
#' p1<-10
#' p2<-100
#' p3<-300
#' ## generation of covariance matrices and feature matrices
#' S1<-riwish(p1,diag(p1))
#' S2<-riwish(p2,diag(p2))
#' S3<-riwish(p3,diag(p3))
#' X1<-matrix(rmvnorm(n*p1,rep(0,p1),S1),n,p1)
#' X2<-matrix(rmvnorm(n*p2,rep(0,p2),S2),n,p2)
#' X3<-matrix(rmvnorm(n*p3,rep(0,p3),S3),n,p3)
#' ## sparsity and generation of betas
#' s2<-p2*0.3
#' s3<-p3*0.01
#' non.zero2<-sample(1:p2,s2)
#' non.zero3<-sample(1:p3,s3)
#' b1<-rnorm(10,0,2.5)
#' b2<-rep(0,p2)
#' b2[non.zero2]<-rnorm(s2)
#' b3<-rep(0,p3)
#' b3[non.zero3]<-rnorm(s3)
## generation of responce
#' mu<-X1%*%b1+X2%*%b2+X3%*%b3
#' y<-rnorm(n,mu,sd=0.5)
#' ## standardize
#' y<-scale(y)
#' X1<-scale(X1)
#' X2<-scale(X2)
#' X3<-scale(X3)
## calculation of gram matrices
#' G1 <- X1 %*% t(X1); G2 <- X2 %*% t(X2); G3 <- X3 %*% t(X3)
#' ## make lists
#' G <- list(G1, G2, G3)
#' X <- list(X1, X2, X3)
#' 
#' ### RUN SBR/SSBR ###
#' 
#' # 1) SBR with the ML lambda-estimator
#' 
#' model1 <- sbr(y = y, X = X,G = G, estimator = 'ML') 
#' 
#' # 2) relaxed SSBR with the ML lambda-estimator using block-matrix computations for the variances of X3 (since p3=300)
#' 
#' model2 <- sbr(y = y, X = X,G = G, estimator = 'ML', sparsify = TRUE, p.threshold = 100, cov.blocks = 100) 
#' 
#' # 3) SSBR with the ML lambda-estimator using block-matrix computations for the variances of X3 (since p3=300) and control equal to log(n) for the effect of sample size
#' 
#' model3 <- sbr(y = y, X = X, G = G, estimator = 'ML', sparsify = TRUE, relaxed = FALSE, p.threshold = 100, cov.blocks = 100, sparse.control = log(n))
#' 
#' # 4) parallel computing for the configuration of model3
#' 
#' cores <- detectCores() - 1
#' cores <- makeCluster(cores)
#' registerDoParallel(cores)
#' model4 <- sbr(y = y, X = X, G = G, estimator = 'ML', parallel = TRUE, cl = cores, sparsify = TRUE, p.threshold = 100, cov.blocks = 50, sparse.control = log(n))
#' stopCluster(cores)
#'
#' ### EXTRACTING OUTPUT FROM A MODEL ###
#'
#' coef(model3)    # SBR/SSBR coefficients (or alternatively model3$coeffients)
#' model3$lambda   # vector of lambdas  
#' model3$sigma2   # error variance  
#' model3$duration # runtime

sbr <- 
function (y, X, trX, G, estimator = "MAP", sparsify = FALSE, relaxed = TRUE, sparse.control = 1, 
    p.threshold = 5000, cov.blocks = 1000, parallel = FALSE, 
    cl, L.optim = 10^-4, U.optim = 10^4) 
{
    start.time <- proc.time()
    if (is.matrix(X) == FALSE & is.list(X) == FALSE) {
        stop("X must be either a matrix (one data-source) or a list (multiple data-sources)")
    }
    if (is.matrix(G) == FALSE & is.list(G) == FALSE) {
        stop("G must be either a matrix (one data-source) or a list (multiple data-sources)")
    }
    if (is.matrix(G) == TRUE) {
        K <- 1
    }
    if (is.list(G) == TRUE) {
        K <- length(G)
    }
    n <- length(y)
    I <- diag(n)
    ty <- t(y)
    if (missing(trX)) {
        if (K == 1) {
            tX <- t(X)
        }
        if (K != 1) {
            tX <- lapply(X, t)
        }
    }
    else {
        tX <- trX
    }
    optim.lambda <- function(args, y, G, estimator, lambda.MAP.prior) {
        if (K == 1) {
            lambda <- args
            lambdaG <- lambda^(-1) * G
        }
        else {
            lambda <- c()
            lambdaG <- list()
            for (j in 1:K) {
                lambda[j] <- args[j]
                lambdaG[[j]] <- lambda[j]^(-1) * G[[j]]
            }
            lambdaG <- Reduce("+", lambdaG)
        }
        mat <- I + lambdaG
        inv.mat <- chol2inv(chol(mat))
        if (estimator == "CV") {
          opt.lambda <- log(ty %*% inv.mat %*% diag(diag(inv.mat^(-2))) %*% 
                              inv.mat %*% y)
        }
        if (estimator == "ML") {
          opt.lambda <- 0.5 * determinant(mat, logarithm = TRUE)$modulus + 
            n/2 * log(0.5 * ty %*% inv.mat %*% y)
        }
        if (estimator == "MAP") {
          opt.lambda <- 0.5 * determinant(mat, logarithm = TRUE)$modulus + 
            n/2 * log(0.5 * ty %*% inv.mat %*% y) + sum(lambda/lambda.MAP.prior)
        }
        opt.lambda
    }
    if (estimator == "CV") {
        lambda <- optim(rep(1,K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "CV")$par
    }
    if (estimator == "ML") {
        lambda <- optim(rep(1,K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "ML")$par
    }
    if (estimator == "MAP") {
        lambda.star <- optim(rep(1,K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "CV")$par
        lambda <- optim(rep(1,K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "MAP", lambda.MAP.prior = lambda.star)$par
    }
    inv.lambda <- lambda^(-1)
    if (K == 1) {
        GL <- inv.lambda * G
    }
    else {
        Glambda <- list()
        for (j in 1:K) {
            Glambda[[j]] <- inv.lambda[j] * G[[j]]
        }
        GL <- Reduce("+", Glambda)
    }
    IGL <- I + GL
    invIGL <- chol2inv(chol(IGL))
    yGL <- y - invIGL %*% GL %*% y
    if (K == 1) {
        beta <- inv.lambda * tX %*% yGL
    }
    else {
        beta <- list()
        for (j in 1:K) {
            beta[[j]] <- inv.lambda[j] * tX[[j]] %*% yGL
        }
        beta <- Reduce(c, beta)
    }
    RSS <- ty %*% invIGL %*% y
    sigma2 <- RSS/(n + 2)
    constant <- c(sqrt(n/RSS))
    message("\n", paste("SBR done."))
    if (sparsify == TRUE) {
        if (K == 1) {
            X <- list(X)
            tX <- list(tX)
        }
        source.weight <- lambda/sum(lambda)
        abs.beta <- abs(beta)
        p.source <- as.numeric(lapply(X, ncol))
        kappa <- (1/abs.beta)^rep(source.weight, times = p.source)
        if(relaxed == FALSE){
          E <- eigen(invIGL)
          sqrt.invIGL <- E$vectors %*% sqrt(diag(E$values)) %*% t(E$vectors)
          sqrt.inv.lambda <- sqrt(inv.lambda)
          sqrt.lambda <- sqrt(lambda)
          M<-list()
          for(j in 1:K){
            M[[j]] <- sqrt.invIGL %*% X[[j]] * sqrt.inv.lambda[j]
          }
          M <- Reduce(cbind, M)
          E <- fast.svd(M)
          options(warn = -1)
          hatD <- sqrt(E$d^2/(1-E$d^2))
          options(warn = 0)
          tV <- t(E$v)
          discard <- which(is.nan(hatD))
          if(length(discard) != 0) {
            hatD <- hatD[-discard]
            tV <- tV[-discard, ]
            l.discard <- length(discard)
            warning(paste("Discarded", l.discard, "negative singular values.", 
                          sep = " "), call. = FALSE)
          }
          subsets <- matrix(NA, K, 2)
          subsets[, 1] <- c(1, 1 + cumsum(p.source)[-K])
          subsets[, 2] <- cumsum(p.source)
          DVL <- list()
          y.part1 <- list()
          y.part2 <- list()
          p <- sum(p.source)
          A <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
          diagA <- c(1/(kappa))
          diag(A) <- diagA
          for(j in 1:K){
            sub <- subsets[j, 1]:subsets[j, 2]
            DVL[[j]] <- diag(hatD) %*% tV[, sub] * sqrt.lambda[j]
            y.part1[[j]] <- DVL[[j]] %*% beta[sub] 
            y.part2[[j]] <- sqrt.lambda[j] * beta[sub]
          }
          DVL <- Reduce(cbind, DVL)
          X.part1 <- Matrix(0, nrow = nrow(DVL), ncol = p)
          X.part1 <- DVL %*% A
          X.part1 <- as.matrix(X.part1)
          X.part2 <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
          diag(X.part2) <- rep(sqrt.lambda, times = p.source)*diagA
          y.star <- constant * c(Reduce('+', y.part1), Reduce('c', y.part2))
          n.star <- length(y.star)
          y.mean <- mean(y.star)
          y.sd <- sd(y.star)
          y.scale <- (y.star - y.mean)/y.sd
          X.star <- constant * rBind(X.part1, X.part2)
          sparse.beta <- diagA * as.numeric(coef(glmnet(X.star, y.scale, lambda = 1/n.star, standardize = TRUE, intercept = FALSE)))[-1]  * y.sd
        } 
          if(relaxed == TRUE) {
              if (parallel == FALSE) {
                ind.small <- which(p.source <= p.threshold)
                ind.big <- which(p.source > p.threshold)
                var.beta <- list()
                length(var.beta) <- K
                D <- diag(1, cov.blocks)
                for (j in ind.small) {
                  var.beta[[j]] <- diag(inv.lambda[j] * (diag(p.source[j]) - 
                                                           tX[[j]] %*% invIGL %*% X[[j]] * inv.lambda[j]))
                }
                for (j in ind.big) {
                  col1 <- seq(1, (p.source[j] - cov.blocks), by = cov.blocks)
                  col2 <- seq(cov.blocks, p.source[j], by = cov.blocks)
                  d <- length(col2)
                  if (length(col1) != d) {
                    col1[d] <- (p.source[j] - cov.blocks) + 1
                  }
                  col2[d] <- p.source[j]
                  for (k in 1:(d - 1)) {
                    flush.console()
                    R <- col1[k]:col2[k]
                    var.beta[[j]][R] <- diag(inv.lambda[j] * (D - 
                                                                tX[[j]][R, ] %*% invIGL %*% X[[j]][, R] * 
                                                                inv.lambda[j]))
                    cat("\r", paste("Source", j, "sigma block-computation:", 
                                    sep = " "), paste(round(k/(d - 1) * 100), 
                                                      "%", sep = ""))
                  }
                  k <- d
                  R <- col1[k]:col2[k]
                  var.beta[[j]][R] <- diag(inv.lambda[j] * (diag(1, 
                                                                 col2[k] - col2[k - 1]) - tX[[j]][R, ] %*% invIGL %*% 
                                                              X[[j]][, R] * inv.lambda[j]))
                  cat("\n", sep = " ")
                }
                var.beta <- Reduce(c, var.beta)
              }
              if (parallel == TRUE) {
                var.function1 <- function(X, inverse.lambda, A) {
                  tX <- t(X)
                  Id <- diag(ncol(X))
                  variance <- diag(inverse.lambda * (Id - tX %*% 
                                                       A %*% X * inverse.lambda))
                  variance
                }
                var.function2 <- function(X, inverse.lambda, A, n.blocks) {
                  tX <- t(X)
                  n.p <- ncol(X)
                  variance <- c()
                  if ((n.p - n.blocks) > n.blocks) {
                    Id <- diag(n.blocks)
                    col1 <- seq(1, (n.p - n.blocks), by = n.blocks)
                    col2 <- seq(n.blocks, n.p, by = n.blocks)
                    d <- length(col2)
                    if (length(col1) != d) {
                      col1[d] <- (n.p - n.blocks) + 1
                    }
                    col2[d] <- n.p
                    for (k in 1:(d - 1)) {
                      R <- col1[k]:col2[k]
                      variance[R] <- diag(inverse.lambda * (Id - 
                                                              tX[R, ] %*% A %*% X[, R] * inverse.lambda))
                    }
                    k <- d
                    R <- col1[k]:col2[k]
                    variance[R] <- diag(inverse.lambda * (diag(1, 
                                                               col2[k] - col2[k - 1]) - tX[R, ] %*% A %*% 
                                                            X[, R] * inverse.lambda))
                  }
                  else {
                    col1 <- c(1, n.blocks + 1)
                    col2 <- c(n.blocks, n.p)
                    for (k in 1:2) {
                      R <- col1[k]:col2[k]
                      variance[R] <- diag(inverse.lambda * (diag(1, 
                                                                 length(R)) - tX[R, ] %*% A %*% X[, R] * 
                                                              inverse.lambda))
                    }
                  }
                  variance
                }
                environment(var.function1) <- .GlobalEnv
                environment(var.function2) <- .GlobalEnv
                var.beta <- list()
                n.cl <- length(cl)
                Xbatches <- list()
                n.col <- rep(NA, K)
                for (j in 1:K) {
                  if (p.source[j] > n.cl) {
                    d <- 1:p.source[j]
                    batches <- split(d, ceiling(seq_along(d)/(p.source[j]/n.cl)))
                    Xbatches[[j]] <- list()
                    for (i in 1:n.cl) {
                      Xbatches[[j]][[i]] <- X[[j]][, batches[[i]]]
                    }
                    n.col[j] <- ncol(Xbatches[[j]][[1]])
                    if (n.col[j] <= p.threshold) {
                      assign(paste("Var.beta"), Reduce(c, parLapply(cl, 
                                                                    as.array(Xbatches[[j]]), var.function1, 
                                                                    inverse.lambda = inv.lambda[j], A = invIGL)))
                    }
                    else {
                      assign(paste("Var.beta"), Reduce(c, parLapply(cl, 
                                                                    as.array(Xbatches[[j]]), var.function2, 
                                                                    inverse.lambda = inv.lambda[j], A = invIGL, 
                                                                    n.blocks = cov.blocks)))
                    }
                  }
                  else {
                    Var.beta <- var.function1(X[[j]], inverse.lambda = inv.lambda[j], 
                                              A = invIGL)
                  }
                  var.beta[[j]] <- Var.beta
                }
                var.beta <- Reduce(c, var.beta)
              }
            sparse.threshold <- 1/constant^2 * var.beta * kappa * sparse.control
            sparse.ind <- which(abs.beta > sparse.threshold)
            sparse.beta <- rep(0, length(beta))
            sparse.beta[sparse.ind] <- beta[sparse.ind] - (beta[sparse.ind]/abs.beta[sparse.ind]) * 
              sparse.threshold[sparse.ind]
          }
        message("\n", paste("SSBR done."))
    }
    diff.time <- proc.time() - start.time
    if (sparsify == FALSE) {
        beta <- as.matrix(beta)
        colnames(beta) <- 'beta'
        Results <- list(coefficients = beta, sigma2 = sigma2, 
                          lambda = lambda, lambdaEstimator = estimator, duration = diff.time)
        } else {
        beta <- cbind(beta, Matrix(sparse.beta))
        colnames(beta) <- c('beta', 'sparse.beta')
        Results <- list(coefficients = beta, sigma2 = sigma2, 
                          lambda = lambda, lambdaEstimator = estimator, duration = diff.time)
        }
    class(Results) <- "sbr"
    Results
}
