# 5 February 2018

## Functions defined below will be used to compute the covariance matrix and the basis functions



## Matern kernel with smoothness parameter 3/2 (from constrKriging)##
Kernel <- function(x, xp, sig, theta)
{
  (sig^2)*(1+(sqrt(3)*abs(x-xp)/theta))*exp(-sqrt(3)*abs(x-xp)/theta)
}



# derivative with respect to the first variable
kp1 <- function(x, xp, sig, theta)
{
  sig^2*(sqrt(3)/theta*sign(x-xp)*exp(-sqrt(3)*(abs(x-xp)/theta))*(-sqrt(3)/theta*abs(x-xp)))
}



# derivative with respect to the second variable
kp2 <- function(x, xp, sig, theta)
{
  -kp1(x,xp, sig, theta)
}



# second derivative with respect to the first and second variables
kpp <- function(x, xp, sig, theta)
  {
    sig^2*((3/theta^2)*exp(-sqrt(3)/theta*abs(x-xp))*(1-sqrt(3)/theta*abs(x-xp)))
  }


# # Gaussian covariance kernel
# Kernel <- function(x, xp, sig, theta){
#   (sig^2)*exp(-(x-xp)^2/(2*theta^2))
# }
# # derivative of the Gaussian kernel with respect to the first variable
# kp1 <- function(x, xp, sig, theta){
#   -(x-xp)/(theta^2)*Kernel(x,xp, sig, theta)
# }
# # derivative with respect to the second variable
# kp2 <- function(x, xp, sig, theta){
#   -kp1(x,xp, sig, theta)
# }
# # second derivative with respect to the 1er and second variables
# kpp <- function(x, xp, sig, theta){
#   (1/(theta^2))*Kernel(x, xp, sig, theta)*(1-(x-xp)^2/(theta^2))
# }



## Covariance matrix ##
covariance_matrix=function(my_knots,sig=1,theta=1)
{
  
  N=length(my_knots)
  CM <- matrix(0, nrow = N + 1, ncol = N + 1)
  
  for(i in 2:(N + 1))
  {
    
    for(j in 2:(N + 1))
    {
       CM[i,j]=kpp(my_knots[i - 1], my_knots[j - 1],sig,theta)
     }
   }
  
  CM[1, 1] <- Kernel(0, 0, sig, theta)
  for(j in 2:(N + 1))
  {
    CM[1, j] <- kp1(0, my_knots[j - 1], sig, theta)
    CM[j, 1] <- kp2(my_knots[j - 1], 0, sig, theta)
  }
    
  return(CM)
 }

# ### Basis function for monotone functions ###
# kappa.am <- function(x, my_knots)
# {
# 
#   N=length(my_knots)
#   k=rep(0,N)
#   p <- rep(0, N)
# 
#   h <- function(x)
#   {
#     if(abs(x) > 1)
#     {
#       return(0)
#     } else {
#       return(1 - abs(x))
#     }
#   }
# 
#   for(j in 1:N)
#   {
#     k[j] <- h((N - 1) * (x - my_knots[j]))
#     ph <- function(t)
#     {
#       return(h((N - 1) * (t - my_knots[j])))
#     }
#     phv <- Vectorize(ph, "t")
#     p[j] <- integrate(phv, lower = 0, upper = x, rel.tol = 1e-10)$value
#     # Set the tolerance at this to avoid the error the integral is probably divergent
#     # and the error roundoff error was detected
#     # and the error extremely bad integrand behaviour
#   }
# 
#   return(p)
# 
# }
# 
# ### Device the basis matrix ###
# basis.am <- function(x, my_knots)
# {
#   phi = matrix(0, length(x), length(my_knots))
#   
#   for(i in 1:length(x))
#   {
#     phi[i, ] = kappa.am(x[i], my_knots)
#   }
#   
#   return(cbind(1, phi))
# }

  
  
### Basis function for monotone functions ###
kappa=function(x,my_knots)
  
{
  
  N=length(my_knots)
  int_length <- (d - c)/(N - 1)
  
  k=rep(0,N)
  
  i=max(which(my_knots<=x))
  
  if(i==1)
    
  {
    
    k[1]=x-0.5*(x^2)/int_length
    
    k[2]=x-my_knots[2]*x/int_length+0.5*x^2/int_length
    
  }
  
  if(i==2)
    
  {
    
    k[1]=int_length/2
    
    k[2]=int_length/2+(x-my_knots[2])*(1+my_knots[2]/int_length)-0.5*(x^2-my_knots[2]^2)/int_length
    
    k[3]=(x-my_knots[2])*(1-my_knots[3]/int_length)+0.5*(x^2-my_knots[2]^2)/int_length
    
  }
  
  if(i==N-1)
    
  {
    
    k[1]=int_length/2
    
    k[2:(N-2)]=int_length
    
    k[N-1]=int_length/2+(x-my_knots[N-1])*(1+my_knots[N-1]/int_length)-0.5*(x^2-my_knots[N-1]^2)/int_length
    
    k[N]=(x-my_knots[N-1])*(1-my_knots[N]/int_length)+0.5*(x^2-my_knots[N-1]^2)/int_length
    
  }
  
  if(i==N)
    
  {
    
    k[1]=int_length/2
    
    k[2:(N-1)]=int_length
    
    k[N]=int_length/2
    
  }
  
  if(i!=1 && i!=2 && i!=N-1 && i!=N)
    
  {
    
    k[1]=int_length/2
    
    k[2:(i-1)]=int_length
    
    k[i]=int_length/2+(x-my_knots[i])*(1+my_knots[i]/int_length)-0.5*(x^2-my_knots[i]^2)/int_length
    
    k[i+1]=(x-my_knots[i])*(1-my_knots[i+1]/int_length)+0.5*(x^2-my_knots[i]^2)/int_length
    
  }
  
  return(k)
  
}



### Device the basis matrix ###
basis <- function(x, my_knots)
{
  phi = matrix(0, length(x), length(my_knots))
  
  for(i in 1:length(x))
  {
    phi[i, ] = kappa(x[i], my_knots)
  }
  
  return(cbind(1, phi))
}


# c <- 0
# d <- 1
# my_knots=seq(0,1, length.out = 5)
# x=seq(0,1,by=0.05)
# phi <- matrix(0, length(x), length(my_knots))
# for(i in 1:length(x))
# {
#   phi[i, ]=kappa(x[i], my_knots)
# }
# plot(x, phi[, 1],type="l",ylim=c(0,0.25),lwd=2)
# for(i in 2:length(my_knots))
# {
#   lines(x, phi[, i],col=i,lwd=2)
# }
# plot(x, phi[, 2],type="l",ylim=c(0, 1),lwd=2)


HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = as.vector(q + epsilon * p)
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # print(c(current_U, proposed_U, exp(current_U-proposed_U), current_K, proposed_K))
  # print(summary(current_p))
  # print(summary(p))
  # print(summary(current_q))
  # print(summary(q))
  
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q)  # accept
  }
  else
  {
    return (current_q)  # reject
  }
}


# ### derivative of basis function for monotone functions ###
happa <- function(x, my_knots)
{
  
  N=length(my_knots)
  k=rep(0,N)
  p <- rep(0, N)
  
  h <- function(x)
  {
    if(abs(x) > 1)
    {
      return(0)
    } else {
      return(1 - abs(x))
    }
  }
  
  for(j in 1:N)
  {
    k[j] <- h(((N - 1)/(d - c)) * (x - my_knots[j]))
  }
  
  return(k)
  
}

### Device the hasis matrix ###
hasis <- function(x, my_knots)
{
  hi = matrix(0, length(x), length(my_knots))
  
  for(i in 1:length(x))
  {
    hi[i, ] = happa(x[i], my_knots)
  }
  
  return(cbind(0, hi))
}


## fast matrix inversion
strassenInv <- function(A){
  
  div4 <- function(A, r){
    A <- list(A)
    A11 <- A[[1]][1:(r/2),1:(r/2)]
    A12 <- A[[1]][1:(r/2),(r/2+1):r]
    A21 <- A[[1]][(r/2+1):r,1:(r/2)]
    A22 <- A[[1]][(r/2+1):r,(r/2+1):r]
    A <- list(X11=A11, X12=A12, X21=A21, X22=A22)
    return(A)
  }
  
  if (nrow(A) != ncol(A)) 
  { stop("only square matrices can be inverted") }
  
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  if ( (is.wholenumber(log(nrow(A), 2)) != TRUE) || (is.wholenumber(log(ncol(A), 2)) != TRUE) )
  { stop("only square matrices of dimension 2^k * 2^k can be inverted with Strassen method") }
  
  A <- div4(A, dim(A)[1])
  
  R1 <- solve(A$X11)
  R2 <- A$X21 %*% R1
  R3 <- R1 %*% A$X12
  R4 <- A$X21 %*% R3
  R5 <- R4 - A$X22
  R6 <- solve(R5)
  C12 <- R3 %*% R6
  C21 <- R6 %*% R2
  R7 <- R3 %*% C21
  C11 <- R1 - R7
  C22 <- -R6
  
  C <- rbind(cbind(C11,C12), cbind(C21,C22))
  
  return(C)
}


condMVN.am <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE) 
{
  if (missing(dependent.ind)) 
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given)) 
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind, 
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind)) 
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma)) 
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-18))   # change this to 10^-18 from 10^-8 not to get error
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}
