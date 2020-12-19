
#' @title chisquqre-estimating interval of population variance.
#' @description estimating interval of population variance using chi-square statistics.
#' @param x the sample
#' @param alpha the significance level
#' @return interval estimation of population variance
#' @importFrom stats qchisq
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100, 1, 2)
#' chisq.var.test(x, 0.05)
#' }
#' @export
chisq.var.test <- function(x, alpha){
  options(digits = 4)
  result <- list()
  n <- length(x)
  v <- var(x)
  result$conf.int.var <- c(
    (n-1)*v/qchisq(alpha/2, df = n-1, lower.tail = F),
    (n-1)*v/qchisq(alpha/2, df = n-1, lower.tail = T)
  )
  result$conf.int.se <- sqrt(result$conf.int.var)
  result
}

#' @title quantile.
#' @description calculate max, min, median, upper and lower quartile.
#' @param x the sample
#' @return max, min, median, upper and lower quartile
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100, 1, 2)
#' my.fivenum(x)
#' }
#' @export
cal.quant <- function(x){
  x <- sort(x)
  n <- length(x)
  n4 <- floor((n+3)/2)/2
  d <- c(1,n4,(n+1)/2,n+1-n4,n)
  return(0.5*(x[floor(d)]+x[ceiling(d)]))
}
