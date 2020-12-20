#' @title Rayleigh sampler
#' @description Rayleigh sampler using acceptance-rejection method
#' @param n size of the random number
#' @param sigma parameter of rayleigh density
#' @return Rayleigh random number and iteration times
#' @examples
#' \dontrun{
#' n <- 1000
#' sigma <- 1
#' x <- Rayleigh(n, sigma)
#' }
#' @export
RayleighR <- function(n, sigma){
  k <- 0
  j <- 0
  result <- list()
  y <- numeric(n)
  while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- rchisq(1,4)
    if (exp(-x^2/(2*sigma^2)+x/2-sigma^2/8) > u) {
      k <- k + 1
      y[k] <- x
    }
  }
  result$rd <- y
  result$it <- j
  result
}
