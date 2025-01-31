# Some constants:
s1  <- 0.6
s2  <- 0.5
mu1 <- 1.5
mu2 <- 1.2

x <- c(0,0,0,0.1,0.1,0.3,0.3,0.9,0.9,0.9)
y <- c(0,0,1,0,1,1,1,0,1,1)

# The function g to be maximised, partial derivatives dg1, dg2, and gradient
g <- function(beta)
{
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  p <- 1 / (1 + exp(-beta0 + beta1 * x))  # logistic regression
  
  log_likelihood <- sum(y * log(p) + (1 - y) * log(1 - p))
  
  return(log_likelihood)
}


gradient <- function(beta)
{
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  p <- 1 / (1 + exp(-(beta0 + beta1 * x)))  # Logistic function
  
  grad_beta0 <- sum(y - p)
  grad_beta1 <- sum((y - p) * x)
  
  return(c(grad_beta0, grad_beta1))
}


# Produce a contour plot; define first a grid where function is evaluated
x1grid <- seq(-2, 2.5, by=0.05)
x2grid <- seq(-2, 3, by=0.05)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- g(c(x1grid[i], x2grid[j]))
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
contour(x1grid, x2grid, mgx, nlevels=34)  # Note: For other functions g, you might need to choose another nlevels-value to get a good impression of the function 

#Steepest ascent function:
steepestasc <- function(x0, eps=1e-8, alpha0=1)
{
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  while(conv>eps)
  {
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*gradient(xt1)
    while (g(xt)<g(xt1))
    {
      alpha <- alpha/2
      xt    <- xt1 + alpha*gradient(xt1)
    }
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- sum((xt-xt1)*(xt-xt1))
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  xt
}
steepestasc(c(-1.5, 2))
