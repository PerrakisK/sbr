#' Function gram
#' 
#' Function for calculating the (inner-product) Gram matrix that allows for block-matrix multiplication.
#' @author Konstanstinos Perrakis \email{konstantinos.perrakis@dzne.de}
#' @author Sach Mukherjee \email{sach.mukherjee@dzne.de}
#' @references Perrakis, K., Mukherjee, S and the Alzheimers Disease Neuroimaging Initiative. (2018) Scalable Bayesian regression in high dimensions with multiple data sources. \url{https://arxiv.org/pdf/1710.00596.pdf}
#' @seealso \code{\link{gram.parallel}}
#' @param X a standardized feature matrix.
#' @param trX (optional) the transpose matrix of X.
#' @param block logical, block matrix computation is performed when TRUE, default option is FALSE.
#' @param block.size used when block = TRUE. Default option is 1000, i.e. blocks of dimensionality 1000 x 1000.
#' @param show.progress logical, when TRUE (and block = TRUE) the progress of the calculations is reported on the console.
#' @return Returns the inner-product Gram matrix.
#' @export
#' @examples 
#' X <- matrix(rnorm(100 * 300), 100, 300)
#' G0 <- gram(X)                             # usual matrix multiplication                          
#' G1 <- gram(X, block=TRUE,block.size=100)  # block matrix multiplication
#' all.equal(G0, G1)  

gram <-
function (X, trX, block = FALSE, block.size = 1000, show.progress = FALSE) 
{
    if (block == FALSE) {
        if (missing(trX)) {
            tX <- t(X)
        }
        else {
            tX <- trX
        }
        G <- X %*% tX
    }
    if (block == TRUE) {
        n <- dim(X)[1]
        p <- dim(X)[2]
        if (block.size == 1) {
            stop("block.size must be greater than 1")
        }
        if (block.size >= p) {
            stop("block.size must me smaller than the number of columns of X")
        }
        if (missing(trX)) {
            tX <- t(X)
        }
        else {
            tX <- trX
        }
        col1 <- seq(1, (p - block.size), by = block.size)
        col2 <- seq(block.size, p, by = block.size)
        d <- length(col2)
        if (length(col1) != d) {
            col1[d] <- (p - block.size) + 1
        }
        col2[d] <- p
        G <- list()
        for (k in 1:d) {
            G[[k]] <- matrix(NA, n, n)
            R <- col1[k]:col2[k]
            G[[k]] <- X[, R] %*% tX[R, ]
            flush.console()
            if(show.progress == TRUE) cat("\r", paste("Gram matrix block-computation:", 
                sep = " "), paste(round(k/d * 100), "%", sep = ""))
        }
        G <- Reduce("+", G)
    }
    G
}
