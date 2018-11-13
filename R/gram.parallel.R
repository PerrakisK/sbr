#' Function gram.parallel
#'
#' Function for calculating the (inner-product) Gram matrix that allows for block-matrix multiplication performed in parallel.
#' @author Konstanstinos Perrakis \email{konstantinos.perrakis@dzne.de}
#' @author Sach Mukherjee \email{sach.mukherjee@dzne.de}
#' @references Perrakis, K., Mukherjee, S. and the Alzheimers Disease Neuroimaging Initiative. (2018) Scalable Bayesian regression in high dimensions with multiple data sources. \url{https://arxiv.org/pdf/1710.00596.pdf}
#' @param X a standardized feature matrix.
#' @param cl the number of cores to use. Must be provided by the user. 
#' @param ... additional arguments passed from function \code{\link{gram}}. One can additionanly use block matrix multiplication within each core. Warning: in this case argument block.size must not exceed floor(ncol(X)/number of cores) - 1.
#' @return Returns the inner-product Gram matrix.
#' @seealso \code{\link{gram}}
#' @export
#' @examples
#' X <- matrix(rnorm(100 * 300), 100, 300)
#' G0 <- gram(X)                        # usual matrix multiplication
#' cores <- detectCores() - 1
#' cores <- makeCluster(cores)
#' registerDoParallel(cores)
#' G1 <- gram.parallel(X, cl = cores)   # parallel matrix multiplication
#' stopCluster(cores)
#' all.equal(G0, G1)

gram.parallel <-
  function (X, cl, ...)
  {
    p <- ncol(X)
    d <- 1:p
    n.cl <- length(cl)
    batches <- split(d, ceiling(seq_along(d) / (p / n.cl)))
    Xbatches <- list()
    for (j in 1:n.cl) {
      Xbatches[[j]] <- X[, batches[[j]]]
    }
    G <-
      Reduce("+", parLapply(cl, Xbatches, gram, ...))
    G
  }
