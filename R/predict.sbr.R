#' Predict SBR/SSBR
#' 
#' Predict S3 method for objects of class 'sbr'. 
#' @author Konstanstinos Perrakis \email{konstantinos.perrakis@dzne.de}
#' @author Sach Mukherjee \email{sach.mukherjee.dzne.de}
#' @references Perrakis, K. and Mukherjee, S. (2018) Scalable Bayesian regression in high dimensions with multiple data sources. \url{https://arxiv.org/pdf/1710.00596.pdf}
#' @param object an object of class 'sbr'.
#' @param newdata a (standardized) data matrix from which to predict.
#' @param coef choose whether to use the SBR beta estimates (default option "sbr") or the SSBR beta estimates (option "ssbr").
#' @return Returns a vector of predictions.
#' @seealso \code{\link{predict}}
#' @export
#' @examples
#' y <- rnorm(100)
#' X <- matrix(rnorm(100*300), 100, 300)
#' G <- gram(X)
#' model <- sbr(y = y, X = X, G = G, sparsify = TRUE)
#' predict1 <- predict(model, X[1:10, ], coef = 'sbr')
#' predict2 <- predict(model, X[1:10, ], coef = 'ssbr')

predict.sbr <-
function (object, newdata, coef = "sbr") 
{
  beta <- object$coefficients
  if ((coef == "ssbr") & (ncol(beta) == 1)) 
    stop("You have requested SSBR coefficients from an SBR fit!")
  if (coef == "sbr") 
    prediction <- newdata %*% beta[, 1]
  if (coef == "ssbr") {
    prediction <- newdata %*% beta[, 2]
  }
  prediction
}
