rm(list = ls())

library('optimParallel')
library('zoo')
library('survival')

closeAllConnections()

# Set up for parallel optimization
cl = makeCluster(spec = detectCores() - 1, outfile = '')
setDefaultCluster(cl = cl)

N = 200  # Number of subjects
T = 25  # Time interval

# Parameters
a12 = 6; a13 = 3; a21 = 1.5; a23 = 4; a31 = 1; a32 = 2;
beta12 = c(1.25); beta13 = c(1.2); beta21 = c(-2.0); beta23 = c(-4.5); beta31 = c(-1.0); beta32 = c(-0.5)

p = c(1/3, 1/3, 1/3)

g = function(x) {
  y = x
  y[y < 100] = log2(1+2^y[y < 100])
  return(y)
}

# Transition Probability Functions
P_1 = function(s, t, x, a12, a13, a21, a23, a31, a32, beta12, beta13, beta21, beta23, beta31, beta32) {
  xb12 = (x %*% beta12)[1]
  xb13 = (x %*% beta13)[1]
  xb21 = (x %*% beta21)[1]
  xb23 = (x %*% beta23)[1]
  xb31 = (x %*% beta31)[1]
  xb32 = (x %*% beta32)[1]
  
  A11 = 1
  A12 = - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { (1+t^a12)^g(xb12)*(1+t^a13)^g(xb13)*((a21*t^(a21-1))/(1+t^a21)*g(xb21)) + 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ((a21*((s+t)/2)^(a21-1))/(1+((s+t)/2)^a21)*g(xb21)) * ((log(1+((s+t)/2)^a12)-log(1+s^a12))/(log(1+t^a12)-log(1+s^a12))) }
  A13 = - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { (1+t^a12)^g(xb12)*(1+t^a13)^g(xb13)*((a31*t^(a31-1))/(1+t^a31)*g(xb31)) + 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ((a31*((s+t)/2)^(a31-1))/(1+((s+t)/2)^a31)*g(xb31)) * ((log(1+((s+t)/2)^a13)-log(1+s^a13))/(log(1+t^a13)-log(1+s^a13))) }
  b1 = ((1+s^a12)^g(xb12)*(1+s^a13)^g(xb13))/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))
  
  A21 = - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { (1+t^a21)^g(xb21)*(1+t^a23)^g(xb23)*((a12*t^(a12-1))/(1+t^a12)*g(xb12)) }
  A22 = 1 - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ((a12*((s+t)/2)^(a12-1))/(1+((s+t)/2)^a12)*g(xb12)) * -((log(1+((s+t)/2)^a12)-log(1+s^a12))/(log(1+t^a12)-log(1+s^a12))) }
  A23 = - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { (1+t^a21)^g(xb21)*(1+t^a23)^g(xb23)*((a32*t^(a32-1))/(1+t^a32)*g(xb32)) + 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ( ((a12*((s+t)/2)^(a12-1))/(1+((s+t)/2)^a12)*g(xb12)) * -((log(1+((s+t)/2)^a13)-log(1+s^a13))/(log(1+t^a13)-log(1+s^a13))) + ((a32*((s+t)/2)^(a32-1))/(1+((s+t)/2)^a32)*g(xb32))*((log(1+((s+t)/2)^a13)-log(1+s^a13))/(log(1+t^a13)-log(1+s^a13))) ) }
  b2 = 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ((a12*((s+t)/2)^(a12-1))/(1+((s+t)/2)^a12)*g(xb12))*1 + (1+s^a21)^g(xb21)*(1+s^a23)^g(xb23)*((a12*s^(a12-1))/(1+s^a12)*g(xb12)) }
  
  A31 = - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { (1+t^a31)^g(xb31)*(1+t^a32)^g(xb32)*((a13*t^(a13-1))/(1+t^a13)*g(xb13)) }
  A32 = - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { (1+t^a31)^g(xb31)*(1+t^a32)^g(xb32)*((a23*t^(a23-1))/(1+t^a23)*g(xb23)) + 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ( ((a13*((s+t)/2)^(a13-1))/(1+((s+t)/2)^a13)*g(xb13)) * -((log(1+((s+t)/2)^a12)-log(1+s^a12))/(log(1+t^a12)-log(1+s^a12))) + ((a23*((s+t)/2)^(a23-1))/(1+((s+t)/2)^a23)*g(xb23))*((log(1+((s+t)/2)^a12)-log(1+s^a12))/(log(1+t^a12)-log(1+s^a12))) ) }
  A33 = 1 - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ( ((a13*((s+t)/2)^(a13-1))/(1+((s+t)/2)^a13)*g(xb13)) * -((log(1+((s+t)/2)^a13)-log(1+s^a13))/(log(1+t^a13)-log(1+s^a13))) ) }
  b3 = 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ((a13*((s+t)/2)^(a13-1))/(1+((s+t)/2)^a13)*g(xb13))*1 + (1+s^a31)^g(xb31)*(1+s^a32)^g(xb32)*((a13*s^(a13-1))/(1+s^a13)*g(xb13)) }
  
  det_A = A11*A22*A33 - A11*A23*A32 - A12*A21*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31
  p_vec = c( (b3*(A12*A23 - A13*A22))/det_A - (b2*(A12*A33 - A13*A32))/det_A + (b1*(A22*A33 - A23*A32))/det_A,
             (b2*(A11*A33 - A13*A31))/det_A - (b3*(A11*A23 - A13*A21))/det_A - (b1*(A21*A33 - A23*A31))/det_A,
             (b3*(A11*A22 - A12*A21))/det_A - (b2*(A11*A32 - A12*A31))/det_A + (b1*(A21*A32 - A22*A31))/det_A )
  L = length(t)
  P = matrix(p_vec, L, 3)
  P = P / rowSums(P)
  return(P)
}

P_2 = function(s, t, x, a12, a13, a21, a23, a31, a32, beta12, beta13, beta21, beta23, beta31, beta32) {
  xb12 = (x %*% beta12)[1]
  xb13 = (x %*% beta13)[1]
  xb21 = (x %*% beta21)[1]
  xb23 = (x %*% beta23)[1]
  xb31 = (x %*% beta31)[1]
  xb32 = (x %*% beta32)[1]
  
  A11 = 1 - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ((a21*((s+t)/2)^(a21-1))/(1+((s+t)/2)^a21)*g(xb21)) * -((log(1+((s+t)/2)^a21)-log(1+s^a21))/(log(1+t^a21)-log(1+s^a21)))
  A12 = - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { (1+t^a12)^g(xb12)*(1+t^a13)^g(xb13)*((a21*t^(a21-1))/(1+t^a21)*g(xb21)) }
  A13 = - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { (1+t^a12)^g(xb12)*(1+t^a13)^g(xb13)*((a31*t^(a31-1))/(1+t^a31)*g(xb31)) + 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ( ((a21*((s+t)/2)^(a21-1))/(1+((s+t)/2)^a21)*g(xb21))*-((log(1+((s+t)/2)^a23)-log(1+s^a23))/(log(1+t^a23)-log(1+s^a23))) + ((a31*((s+t)/2)^(a31-1))/(1+((s+t)/2)^a31)*g(xb31)) * ((log(1+((s+t)/2)^a23)-log(1+s^a23))/(log(1+t^a23)-log(1+s^a23))) ) }
  b1 = 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ((a21*((s+t)/2)^(a21-1))/(1+((s+t)/2)^a21)*g(xb21))*1 + (1+s^a12)^g(xb12)*(1+s^a13)^g(xb13)*((a21*s^(a21-1))/(1+s^a21)*g(xb21)) }
  
  A21 = - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { (1+t^a21)^g(xb21)*(1+t^a23)^g(xb23)*((a12*t^(a12-1))/(1+t^a12)*g(xb12)) + 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ((a12*((s+t)/2)^(a12-1))/(1+((s+t)/2)^a12)*g(xb12)) * ((log(1+((s+t)/2)^a21)-log(1+s^a21))/(log(1+t^a21)-log(1+s^a21)))}
  A22 = 1
  A23 = - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { (1+t^a21)^g(xb21)*(1+t^a23)^g(xb23)*((a32*t^(a32-1))/(1+t^a32)*g(xb32)) + 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ((a32*((s+t)/2)^(a32-1))/(1+((s+t)/2)^a32)*g(xb32)) * ((log(1+((s+t)/2)^a23)-log(1+s^a23))/(log(1+t^a23)-log(1+s^a23)))}
  b2 = ((1+s^a21)^g(xb21)*(1+s^a23)^g(xb23)) / ((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))
  
  A31 = - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { (1+t^a31)^g(xb31)*(1+t^a32)^g(xb32)*((a13*t^(a13-1))/(1+t^a13)*g(xb13)) + 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ( ((a13*((s+t)/2)^(a13-1))/(1+((s+t)/2)^a13)*g(xb13)) * ((log(1+((s+t)/2)^a21)-log(1+s^a21))/(log(1+t^a21)-log(1+s^a21))) + ((a23*((s+t)/2)^(a23-1))/(1+((s+t)/2)^a23)*g(xb23))*-((log(1+((s+t)/2)^a21)-log(1+s^a21))/(log(1+t^a21)-log(1+s^a21))) ) }
  A32 = - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { (1+t^a31)^g(xb31)*(1+t^a32)^g(xb32)*((a23*t^(a23-1))/(1+t^a23)*g(xb23)) }
  A33 = 1 - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ( ((a23*((s+t)/2)^(a23-1))/(1+((s+t)/2)^a23)*g(xb23)) * -((log(1+((s+t)/2)^a23)-log(1+s^a23))/(log(1+t^a23)-log(1+s^a23))) ) }
  b3 = 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ((a23*((s+t)/2)^(a23-1))/(1+((s+t)/2)^a23)*g(xb23))*1 + (1+s^a31)^g(xb31)*(1+s^a32)^g(xb32)*((a23*s^(a23-1))/(1+s^a23)*g(xb23)) }
  
  det_A = A11*A22*A33 - A11*A23*A32 - A12*A21*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31
  p_vec = c( (b3*(A12*A23 - A13*A22))/det_A - (b2*(A12*A33 - A13*A32))/det_A + (b1*(A22*A33 - A23*A32))/det_A,
             (b2*(A11*A33 - A13*A31))/det_A - (b3*(A11*A23 - A13*A21))/det_A - (b1*(A21*A33 - A23*A31))/det_A,
             (b3*(A11*A22 - A12*A21))/det_A - (b2*(A11*A32 - A12*A31))/det_A + (b1*(A21*A32 - A22*A31))/det_A )
  L = length(t)
  P = matrix(p_vec, L, 3)
  P = P / rowSums(P)
  return(P)
}

P_3 = function(s, t, x, a12, a13, a21, a23, a31, a32, beta12, beta13, beta21, beta23, beta31, beta32) {
  xb12 = (x %*% beta12)[1]
  xb13 = (x %*% beta13)[1]
  xb21 = (x %*% beta21)[1]
  xb23 = (x %*% beta23)[1]
  xb31 = (x %*% beta31)[1]
  xb32 = (x %*% beta32)[1]
  
  A11 = 1 - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ((a31*((s+t)/2)^(a31-1))/(1+((s+t)/2)^a31)*g(xb31)) * -((log(1+((s+t)/2)^a31)-log(1+s^a31))/(log(1+t^a31)-log(1+s^a31))) }
  A12 = - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { (1+t^a12)^g(xb12)*(1+t^a13)^g(xb13)*((a21*t^(a21-1))/(1+t^a21)*g(xb21)) + 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ( ((a21*((s+t)/2)^(a21-1))/(1+((s+t)/2)^a21)*g(xb21))*((log(1+((s+t)/2)^a32)-log(1+s^a32))/(log(1+t^a32)-log(1+s^a32))) + ((a31*((s+t)/2)^(a31-1))/(1+((s+t)/2)^a31)*g(xb31)) * -((log(1+((s+t)/2)^a32)-log(1+s^a32))/(log(1+t^a32)-log(1+s^a32))) ) }
  A13 = - 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { (1+t^a12)^g(xb12)*(1+t^a13)^g(xb13)*((a31*t^(a31-1))/(1+t^a31)*g(xb31))  }
  b1 = 1/((1+t^a12)^g(xb12)*(1+t^a13)^g(xb13))*(t-s)/6 * { 4*(1+((s+t)/2)^a12)^g(xb12)*(1+((s+t)/2)^a13)^g(xb13) * ((a31*((s+t)/2)^(a31-1))/(1+((s+t)/2)^a31)*g(xb31))*1 + (1+s^a12)^g(xb12)*(1+s^a13)^g(xb13)*((a31*s^(a31-1))/(1+s^a31)*g(xb31)) }
  
  A21 = - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { (1+t^a21)^g(xb21)*(1+t^a23)^g(xb23)*((a12*t^(a12-1))/(1+t^a12)*g(xb12)) + 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * (((a12*((s+t)/2)^(a12-1))/(1+((s+t)/2)^a12)*g(xb12)) * ((log(1+((s+t)/2)^a31)-log(1+s^a31))/(log(1+t^a31)-log(1+s^a31))) + ((a32*((s+t)/2)^(a32-1))/(1+((s+t)/2)^a32)*g(xb32))*-((log(1+((s+t)/2)^a31)-log(1+s^a31))/(log(1+t^a31)-log(1+s^a31))))}
  A22 = 1 - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ( ((a32*((s+t)/2)^(a32-1))/(1+((s+t)/2)^a32)*g(xb32)) * -((log(1+((s+t)/2)^a32)-log(1+s^a32))/(log(1+t^a32)-log(1+s^a32))) ) }
  A23 = - 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { (1+t^a21)^g(xb21)*(1+t^a23)^g(xb23)*((a32*t^(a32-1))/(1+t^a32)*g(xb32)) }
  b2 = 1/((1+t^a21)^g(xb21)*(1+t^a23)^g(xb23))*(t-s)/6 * { 4*(1+((s+t)/2)^a21)^g(xb21)*(1+((s+t)/2)^a23)^g(xb23) * ((a32*((s+t)/2)^(a32-1))/(1+((s+t)/2)^a32)*g(xb32))*1 + (1+s^a21)^g(xb21)*(1+s^a23)^g(xb23)*((a32*s^(a32-1))/(1+s^a32)*g(xb32)) }
  
  A31 = - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { (1+t^a31)^g(xb31)*(1+t^a32)^g(xb32)*((a13*t^(a13-1))/(1+t^a13)*g(xb13)) + 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ( ((a13*((s+t)/2)^(a13-1))/(1+((s+t)/2)^a13)*g(xb13)) * ((log(1+((s+t)/2)^a31)-log(1+s^a31))/(log(1+t^a31)-log(1+s^a31))) ) }
  A32 = - 1/((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))*(t-s)/6 * { (1+t^a31)^g(xb31)*(1+t^a32)^g(xb32)*((a23*t^(a23-1))/(1+t^a23)*g(xb23)) + 4*(1+((s+t)/2)^a31)^g(xb31)*(1+((s+t)/2)^a32)^g(xb32) * ( ((a23*((s+t)/2)^(a23-1))/(1+((s+t)/2)^a23)*g(xb23)) * ((log(1+((s+t)/2)^a32)-log(1+s^a32))/(log(1+t^a32)-log(1+s^a32))) )}
  A33 = 1
  b3 = ((1+s^a31)^g(xb31)*(1+s^a32)^g(xb32)) / ((1+t^a31)^g(xb31)*(1+t^a32)^g(xb32))
  
  det_A = A11*A22*A33 - A11*A23*A32 - A12*A21*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31
  p_vec = c( (b3*(A12*A23 - A13*A22))/det_A - (b2*(A12*A33 - A13*A32))/det_A + (b1*(A22*A33 - A23*A32))/det_A,
             (b2*(A11*A33 - A13*A31))/det_A - (b3*(A11*A23 - A13*A21))/det_A - (b1*(A21*A33 - A23*A31))/det_A,
             (b3*(A11*A22 - A12*A21))/det_A - (b2*(A11*A32 - A12*A31))/det_A + (b1*(A21*A32 - A22*A31))/det_A )
  L = length(t)
  P = matrix(p_vec, L, 3)
  P = P / rowSums(P)
  return(P)
}


# Data Simulation
MC = matrix(rep(NA, N * (T + 1)), nrow = N, ncol = (T + 1))
X = matrix(runif(N), ncol=1)
n_covar = ncol(X)

for (i in 1:N) {
  state = sample(c(1,2,3), size=1, prob=p)
  timer = c(0)
  
  while (timer[length(timer)] <= T) {
    tm = timer[length(timer)]
    st = state[length(state)]
    Z1 = rexp(1, rate=1)
    Z2 = rexp(1, rate=1)
    
    if (st == 1) {
      t12 = (exp(Z1 / g((X[i,] %*% beta12)[1]))*(1+tm^a12)-1)^(1/a12)
      t13 = (exp(Z2 / g((X[i,] %*% beta13)[1]))*(1+tm^a13)-1)^(1/a13)
      t_next = min(t12, t13)
      s_next = which.min(c(Inf, t12, t13))
    } else if (st == 2) {
      t21 = (exp(Z1 / g((X[i,] %*% beta21)[1]))*(1+tm^a21)-1)^(1/a21)
      t23 = (exp(Z2 / g((X[i,] %*% beta23)[1]))*(1+tm^a23)-1)^(1/a23)
      t_next = min(t21, t23)
      s_next = which.min(c(t21, Inf, t23))
    } else if (st == 3) {
      t31 = (exp(Z1 / g((X[i,] %*% beta31)[1]))*(1+tm^a31)-1)^(1/a31)
      t32 = (exp(Z2 / g((X[i,] %*% beta32)[1]))*(1+tm^a32)-1)^(1/a32)
      t_next = min(t31, t32)
      s_next = which.min(c(t31, t32, Inf))
    }
    
    timer = c(timer, t_next)
    state = c(state, s_next)
  }
  
  for (j in 1:length(timer)) {
    if (timer[j] > T) break
    MC[i, ceiling(timer[j]) + 1] = state[j]
  }
}

for (i in 1:nrow(MC)) {
  MC[i, ] = na.locf(MC[i, ])
}


# Initial value estimation
A = matrix(NA, N, 6)
for (i in 1:N) {
  Y = MC[i,]
  tt = c(0:T); tt[1] = 0.001
  rate12 = c(); rate13 = c(); rate21 = c(); rate23 = c(); rate31 = c(); rate32 = c();
  st = Y[1]; tm = tt[1]
  
  for (k in 2:length(tt)) {
    st2 = Y[k]; tm2 = tt[k]
    if (st2 != st) {
      rate_val = 1/(log(tm2)-log(tm))
      if (st == 1 && st2 == 2) rate12 = c(rate12, rate_val)
      else if (st == 1 && st2 == 3) rate13 = c(rate13, rate_val)
      else if (st == 2 && st2 == 1) rate21 = c(rate21, rate_val)
      else if (st == 2 && st2 == 3) rate23 = c(rate23, rate_val)
      else if (st == 3 && st2 == 1) rate31 = c(rate31, rate_val)
      else if (st == 3 && st2 == 2) rate32 = c(rate32, rate_val)
      st = st2; tm = tm2
    }
  }
  
  rates = list(rate12, rate13, rate21, rate23, rate31, rate32)
  for (l in 1:6) {
    if (length(rates[[l]]) == 0) A[i, l] = 1e-3
    else A[i, l] = min(mean(rates[[l]]), 7)
  }
}

m1 = colMeans(A, na.rm = TRUE)

long_list = list()
for (i in 1:N) {
  Y = as.numeric(MC[i,])
  tt = 0:T
  covs = as.numeric(X[i, ]) 
  
  for (j in 1:(length(Y) - 1)) {
    long_list[[length(long_list) + 1]] = c(i, tt[j], tt[j+1], Y[j], Y[j+1], covs)
  }
}
long_dat = as.data.frame(do.call(rbind, long_list))
cov_names = paste0('cov', 1:n_covar)
colnames(long_dat) = c('id', 'start', 'stop', 'origin', 'dest', cov_names)
cox_form = as.formula(paste('Surv(start, stop, event) ~', paste(cov_names, collapse=' + ')))

beta_vec <- c()
ordered_trans = list(c(1,2), c(1,3), c(2,1), c(2,3), c(3,1), c(3,2))

for (tr in ordered_trans) {
  dsub <- subset(long_dat, origin == tr[1])
  dsub$event <- as.integer(dsub$dest == tr[2])
  
  fit_cox = tryCatch({ coxph(cox_form, data = dsub) }, error = function(e) NULL)
  
  if (is.null(fit_cox)) {
    beta_vec = c(beta_vec, rep(0, n_covar))
  } else {
    coefs = unname(coef(fit_cox))
    coefs[is.na(coefs)] = 0
    beta_vec = c(beta_vec, coefs)
  }
}

par0 = c(m1, beta_vec)


# Negative log-likelihood function
nLL = function(par, MC, X) {
  sum_ll = 0
  timeseq = 0:T; timeseq[1] = 1e-15
  n_covar = ncol(X)
  
  b12 = par[(7)             : (6 + 1*n_covar)]
  b13 = par[(7 + 1*n_covar) : (6 + 2*n_covar)]
  b21 = par[(7 + 2*n_covar) : (6 + 3*n_covar)]
  b23 = par[(7 + 3*n_covar) : (6 + 4*n_covar)]
  b31 = par[(7 + 4*n_covar) : (6 + 5*n_covar)]
  b32 = par[(7 + 5*n_covar) : (6 + 6*n_covar)]
  
  for (i in 1:nrow(MC)) {
    s_state = MC[i, 1:T]; t_state = MC[i, 2:(T+1)]
    x = as.numeric(X[i, ])
    
    P1_mat = P_1(timeseq[1:T], timeseq[2:(T+1)], x, par[1], par[2], par[3], par[4], par[5], par[6], b12, b13, b21, b23, b31, b32)
    P2_mat = P_2(timeseq[1:T], timeseq[2:(T+1)], x, par[1], par[2], par[3], par[4], par[5], par[6], b12, b13, b21, b23, b31, b32)
    P3_mat = P_3(timeseq[1:T], timeseq[2:(T+1)], x, par[1], par[2], par[3], par[4], par[5], par[6], b12, b13, b21, b23, b31, b32)
    
    P_mat = sweep(P1_mat, MARGIN=1, as.numeric(s_state == 1), `*`) + 
      sweep(P2_mat, MARGIN=1, as.numeric(s_state == 2), `*`) + 
      sweep(P3_mat, MARGIN=1, as.numeric(s_state == 3), `*`)
    
    P_mat[is.na(P_mat) | P_mat <= 0] = 1e-15
    Prod = prod(P_mat[cbind(1:nrow(P_mat), t_state)])
    
    sum_ll = sum_ll - log(Prod)
  }
  
  if (!is.finite(sum_ll) || is.nan(sum_ll)) {sum_ll = 1e100}
  return(sum_ll)
}


clusterExport(cl, c('P_1', 'P_2', 'P_3', 'T', 'g'))

Start_time = Sys.time()

# Primary optimization with L-BFGS-B bounds
fit = try(optimParallel(
  par = par0, 
  fn = nLL, 
  MC = MC, 
  X = X, 
  control = list(maxit=1000), 
  lower = c(rep(1e-4, 6), rep(-Inf, 6*n_covar)), 
  method = 'L-BFGS-B', 
  hessian = TRUE
))

# Fallback to Nelder-Mead if L-BFGS-B fails
if (inherits(fit, 'try-error')) {
  print('L-BFGS-B failed. Falling back to constrOptim (Nelder-Mead)...')
  
  num_params = length(par0)
  num_constraints = 6
  
  UI = matrix(0, nrow = num_constraints, ncol = num_params)
  CI = rep(0, num_constraints)
  
  for (j in 1:6) {
    UI[j, j] = 1
  }
  
  fit = try(constrOptim(
    theta = par0, 
    f = nLL, 
    MC = MC, 
    X = X, 
    ui = UI, 
    ci = CI, 
    control = list(maxit=1000), 
    method = 'Nelder-Mead',
    hessian = TRUE
  ))
}

End_time = Sys.time()
print(difftime(End_time, Start_time, unit='mins'))
print('Simulation run completed.')

stopCluster(cl)


SE = tryCatch({
  sqrt(diag(solve(fit$hessian)))
}, error = function(e) {
  rep(NA, length(fit$par))
})

alpha_labels = paste0('alpha_', c('12', '13', '21', '23', '31', '32'))
beta_labels = c()

for (tr in c('12', '13', '21', '23', '31', '32')) {
  for (i in 1:n_covar) {
    beta_labels = c(beta_labels, paste0('beta_', tr, '_cov', i))
  }
}

results_df = data.frame(
  Parameter = c(alpha_labels, beta_labels),
  Estimate = round(fit$par, 4),
  Std_Error = round(SE, 4)
)

print(results_df)
