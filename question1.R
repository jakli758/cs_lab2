f <- function(x, y) {
  return (sin(x+y) + (x-y)^2 - 1.5*x + 2.5*y + 1)
}

# a. contour plot

# using example code from lecture
x1grid <- -30:80/20
x2grid <- -60:80/20
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx <- dx1*dx2

fx <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    x <- x1grid[i]
    y <- x2grid[j]
    fx[i, j] <- f(x, y)
  }
mfx <- matrix(fx, nrow=dx1, ncol=dx2)


contour(x1grid, x2grid, mfx, nlevels=20, xlab="x", ylab="y", main="Contour Plot for f(x,y)")
contour(x1grid, x2grid, mfx, nlevels=100, xlab="x", ylab="y", main="Contour Plot for f(x,y)")
contour(x1grid, x2grid, mfx, nlevels=200, xlab="x", ylab="y", main="Contour Plot for f(x,y)")

# 3d plot, not asked for in task
# plot is mirrored to better see the maxima (minima in the original data)
mfxcut <- mfx/(mfx < 10)   # set values below -15 to -Inf to avoid them in the plot
persp(x1grid, x2grid, -mfxcut, xlab="x", ylab="y", zlab="z", theta=45, phi=20, zlim=c(-15, 15))

# b. Gradient + Hessian

# get gradient
deriv(expression(sin(x+y) + (x-y)^2 - 1.5*x + 2.5*y + 1), namevec = c("x", "y"))

grad_g <- function(x, y) {
  grad_x <- cos(x + y) + 2 * (x - y) - 1.5
  grad_y <- cos(x + y) - 2 * (x - y) + 2.5
  return(c(grad_x, grad_y))
}

hessian_g <- function(x, y) {
  h11 <- -sin(x+y) + 2
  h12 <- -sin(x+y)
  h21 <- -sin(x+y)
  h22 <- -sin(x+y) + 2
  return(matrix(c(h11, h12, h21, h22), nrow = 2, byrow = TRUE))
}

a <- c(-1,1)
a[1]
hessian_g(a[1], a[1])

# c. custom newton

newton <- function(x0,y0) {
  x_t <- c(x0, y0)
  
  it <- 0
  while(TRUE) {
    x_t1 <- x_t - solve(hessian_g(x_t[1],x_t[2])) %*% grad_g(x_t[1],x_t[2])
    
    # use absolute stopping criterion or a maximum number of iterations
    if (norm(x_t1 - x_t, type="2") < 0.005 | it >= 10000){
      break
    }
    
    x_t <- x_t1
    it <- it + 1
  }
  
  return(c(x_t, it))
}

newton(0,-1) # finds minimum with convergence in 7 iterations, 
# at the same point this is a global minimum if we look at the contour plot


newton(2,-1) # does not find an optimum within 1000 iterations, diverges
newton(3,2) # finds minimum in 7 iterations, but only local
newton(1.5,0.5)
