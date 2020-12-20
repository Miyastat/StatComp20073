#' @title calculate confidence interval.
#' @description A function used to calculate confidence interval.
#' @param n the number of samples.
#' @param alpha significance level.
#' @return confidence border
#' @importFrom stats rlnorm sd qt qchisq rchisq rt cor runif rnorm dnorm cov var qnorm pnorm quantile
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp20073
#' @examples
#' \dontrun{
#' C1 <- replicate(1000, expr = calcCI(n = 20, alpha = .025))
#' C2 <- replicate(1000, expr = calcCI(n = 20, alpha = .975))
#' print(c(exp(mean(C1)), exp(mean(C2)))) ## the confidence interval
#' }
#' @export
calcCI <- function(n, alpha) {
  x <- rlnorm(n, meanlog = 1, sdlog = 2)
  y <- log(x)
  return(sd(y) * qt(alpha, df = n-1) / sqrt(n) + mean(y))
}


#' @title count five method for hypothesis tests.
#' @description count five method for hypothesis tests.
#' @param x a (non-empty) numeric vector of data values.
#' @param y a vector, matrix or data frame with compatible dimensions to x.
#' @return the value 1 (reject null hypothesis) or 0 (do not reject null hypothesis)
#' @examples
#' \dontrun{
#' power <- mean(replicate(1e3, expr={
#'   x <- rnorm(20, 0, 1)
#'   y <- rnorm(20, 0, 1.5)
#'   count5test(x, y)
#' }))
#' }
#' @export
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}

#' @title skewness.
#' @description skewness.
#' @param x a numeric vector, matrix or data frame.
#' @return the skewness of the sample
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' res <- sk(x)
#' }
#' @export
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#' @title Mardia tests of multivariate normality.
#' @description tests of multivariate normality based on multivariate generalizations of skewness and kurtosis proposed by Mardia.
#' @param X the observation.
#' @return the empirical estimates of power
#' @export
Mardia.test <- function(X) {
  n <- nrow(X)
  c <- ncol(X)
  mydata <- X
  for(i in 1:c){
    mydata[,i]<-X[,i]-mean(X[,i])
  }
  sigmah<-t(mydata)%*%mydata/n
  a<-mydata%*%solve(sigmah)%*%t(mydata)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

#' @title A Metropolis sampler using R.
#' @description A Metropolis sampler using R.
#' @param n the number of the samples.
#' @param sigma vector of standard deviations.
#' @param x0 initial value.
#' @param N  the length of the chain.
#' @return a random walk Metropolis sample of size
#' @examples
#' \dontrun{
#' rw <- MetropolisR(4, 2, 25, 2000)
#' plot(rw$x, type = "l")
#' }
#' @export
MetropolisR <- function(n, sigma, x0, N) {
  f <- function(x) {
    return(exp(-abs(x))/2)
  }
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y) / f(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

#' @title The EM algorithm.
#' @description The EM algorithm.
#' @param ep critical value
#' @return The EM estimate
#' @export
est <- function(ep = 1e-7){
  p = 0.2; q = 0.1
  k = 1
  repeat{
    k = k+1
    p[k]=(444*p[k-1]/(2-p[k-1]-2*q[k-1])+507)/2000
    q[k]=(132*q[k-1]/(2-q[k-1]-2*p[k-1])+195)/2000
    if(abs(p[k]-p[k-1])<ep | abs(q[k]-q[k-1])<ep)break
  }
  list(p=p[k], q=q[k], iter=k)
}

#' @title calculate confidence interval.
#' @description A function used to calculate confidence interval.
#' @param x1 A list or atomic vector.
#' @param x2 A list or atomic vector.
#' @param FUN the function to be applied to each element of X.
#' @param FUN.VALUE a (generalized) vector; a template for the return value from FUN. See ‘Details’.
#' @param USE.NAMES logical.
#' @examples
#' \dontrun{
#' set.seed(123)
#' xs <- list(runif(10), runif(20))
#' ws <- list(rpois(10, 2) + 1, rpois(20, 2) + 1)
#' res <- nlapply(xs,ws,weighted.mean, numeric(1))
#' }
#' @export
nlapply <- function(x1,x2, FUN, FUN.VALUE, USE.NAMES = TRUE){
  answer <- Map(FUN, x1,x2)
  vapply(answer, function(x) x, FUN.VALUE = FUN.VALUE)
}

#' @title the Gelman-Rubin method of monitoring convergence of a Metropolis chain.
#' @description the Gelman-Rubin method of monitoring convergence of a Metropolis chain.
#' @param psi  the proposal distribution.
#' @return The Gelman-Rubin statistic
#' @export
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, "var")
  W <- mean(psi.w)
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W
  return(r.hat)
}

#' @title The Metropolis-Hastings sampler.
#' @description The Metropolis-Hastings sampler.
#' @param sigma standard deviation of proposal distribution
#' @param N length of the chain
#' @param X1 mean of proposal distribution
#' @return confidence border
#' @export
normal.chain <- function(sigma, N, X1) {
  ldf <- function(x,u,b) exp(-abs(x-u)/b)/2/b
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rnorm(1, xt, sigma)
    r1 <- ldf(y, 0, 1) * dnorm(xt, y, sigma)
    r2 <- ldf(xt, 0, 1) * dnorm(y, xt, sigma)
    r <- r1 / r2
    if (u[i] <= r) x[i] <- y else
      x[i] <- xt
  }
  return(x)
}

#' @title the maximum number of extreme points.
#' @description counts the maximum number of extreme points of each sample with respect to the range of the other sample.
#' @param x numeric or complex or logical vectors.
#' @param y numeric or complex or logical vectors.
#' @return the maximum number of extreme points
#' @examples
#' \dontrun{
#' x <- rnorm(10, 0, 1)
#' y <- rnorm(10, 0, 1)
#' res <- maxout(x, y)
#' }
#' @export
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

#' @title BCa bootstrap confidence interval.
#' @description BCa bootstrap confidence interval.
#' @param x the observations
#' @param th0 the observed statistic.
#' @param th is the vector of bootstrap replicates.
#' @param stat the function to compute the statistic
#' @param conf the confidence level
#' @return the maximum number of extreme points
#' @export
boot.BCa <-
  function(x, th0, th, stat, conf = .95) {
    x <- as.matrix(x)
    n <- nrow(x)
    N <- 1:n
    alpha <- (1 + c(-conf, conf))/2
    zalpha <- qnorm(alpha)
    z0 <- qnorm(sum(th < th0) / length(th))
    th.jack <- numeric(n)
    for (i in 1:n) {
      J <- N[1:(n-1)]
      th.jack[i] <- stat(x[-i, ], J)
    }
    L <- mean(th.jack) - th.jack
    a <- sum(L^3)/(6 * sum(L^2)^1.5)
    adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
    limits <- quantile(th, adj.alpha, type=6)
    return(list("est"=th0, "BCa"=limits))
  }

#' @title function to compute the statistic.
#' @description function to compute the statistic.
#' @param dat the observed data
#' @param ind the observed data
#' @return a simulated data set of the same form as the observed data
#' @export
theta.boot <- function(dat, ind) {
  x <- dat[ind, 1]
  mean(x)
}
NULL

#' @title the empirical power using  \code{t(v)}
#' @description calculate the empirical power with heavy-tailed t.
#' @param v degrees of freedom
#' @return the empirical power
#' @export
pwr_t = function(v){

  alpha = 0.1
  n = 20
  m = 1e3
  N = length(v)
  pwr = numeric(N)
  cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

  for (j in 1:N) {
    sktests = numeric(m)
    for (i in 1:m) {
      x = rt(n,v[j])
      sktests[i] = as.integer(abs(sk(x))>= cv)
    }
    pwr[j] = mean(sktests)
  }
  se = sqrt(pwr*(1-pwr) / m)
  return(list(pwr = pwr,se = se))
}
