x <- c(0,0,0,0.1,0.1,0.3,0.3,0.9,0.9,0.9)
y <- c(0,0,1,0,1,1,1,0,1,1)

# The function g to be maximised, partial derivatives dg1, dg2, and gradient
g <- function(beta)
{
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  p <- 1 / (1 + exp(-(beta0 + beta1 * x)))  # logistic regression
  
  log_likelihood <- sum(y * log(p) + (1 - y) * log(1 - p))
  
  return(log_likelihood)
}


gradient <- function(beta)
{
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  p <- 1 / (1 + exp(-(beta0 + beta1 * x)))
  
  grad_beta0 <- sum(y - p)
  grad_beta1 <- sum((y - p) * x)
  
  return(c(grad_beta0, grad_beta1))
}


# Produce a contour plot; define first a grid where function is evaluated
x1grid <- seq(-4, 4, by=0.05)
x2grid <- seq(-4, 4, by=0.05)
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
contour(x1grid, x2grid, mgx, nlevels=50)  # Note: For other functions g, you might need to choose another nlevels-value to get a good impression of the function 

#Steepest ascent function:

# TODO: ensure five digits precision
steepestasc <- function(x0, eps=1e-10, alpha0=1, factor=0.5)
{
  func_evals <- 0
  grad_evals <- 0
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  while(conv>eps)
  {
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*gradient(xt1)
    grad_evals <- grad_evals + 1
    
    while (g(xt)<g(xt1))
    {
      func_evals <- func_evals + 2 # for the two evals in while condition
      
      # adjust alpha by given factor
      alpha <- alpha*factor
      
      xt    <- xt1 + alpha*gradient(xt1)
      grad_evals <- grad_evals + 1
    }
    func_evals <- func_evals + 2 # for the two evals in while condition that broke the loop
    
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- sum((xt-xt1)*(xt-xt1))
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  
  cat("Function evaluations:", func_evals, "\n")
  cat("Gradient evaluations:", grad_evals, "\n")
  
  xt
}

# according to gentle:
# "For a function whose contours are ellipses, as the function in Exercise 6.10
# (page 302), for example, the steepest descent steps will zigzag toward the solution" (p. 265)


steepestasc(c(-2, -2))
steepestasc(c(-2, -2), factor=0.4)
steepestasc(c(-2, 2))
steepestasc(c(2.5, 0))
steepestasc(c(-0.2, 1))
steepestasc(c(-0.2, 1), factor=0.4)
contour(x1grid, x2grid, mgx, nlevels=50)


## task c.


optim(c(-0.2, 1), g, gr=gradient, method="BFGS", control=list(fnscale=-1)) # use fnscale to maximize

optim(c(-0.2, 1), g, method="Nelder-Mead",control=list(fnscale=-1))


## d. use glm

model <- glm(y ~ x, family="binomial")
model$iter # the number of iterations of IWLS used.
model$coefficients
