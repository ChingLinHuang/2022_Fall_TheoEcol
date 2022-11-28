library(tidyverse)
library(deSolve)

### Model specification
RM_predation_model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N*(1-(N/K))-(a*N/(1+a*h*N))*P
    dP_dt = e*(a*N/(1+a*h*N))*P-d*P
    return(list(c(dN_dt, dP_dt)))
  })
}

### Model parameters
times <- seq(0, 2000, by = 0.01)
state <- c(N = 5, P = 2)
parms <- c(r = 1.0, K = 7.0, a = 1.3, h = 0.9, e = 0.6, d = 0.5)

### Model application
pop_size <- ode(func = RM_predation_model, times = times, y = state, parms = parms)

### equilibrium
E_np <- with(as.list(parms),
             c(N = d/(a*(e-d*h)),
               P = r/a*(1-d/(a*(e-d*h))/K)*(1+a*h*d/(a*(e-d*h)))))


### Visualize the population dynamics
# population trajectories
plot(c(0, max(times)), c(0, max(pop_size[, c("N", "P")])), type = "n", xlab = "time", ylab = "population size")
lines(N ~ time, data = pop_size, col = "blue") # dynamics of N
lines(P ~ time, data = pop_size, col = "red") # dynamics of P
legend("topright", legend = c("N", "P"), col = c("blue", "red"), lty = 1)

# state-space diagram
max_P <- max(pop_size[ ,"P"])
max_N <- max(pop_size[ ,"N"])
plot(P ~ N, data = pop_size, type = "l", xlim = c(0, max_N), ylim = c(0, max_P*1.1))
points(E_np["P"] ~ E_np["N"], pch = 16) # equilibrium
with(as.list(parms), {
  # ZNGIs of N
  abline(v = 0, col = "blue")
  curve(r/a*(1-x/K)*(1+a*h*x), from = -2, to = K+2, col = "blue", add = T)
  # ZNGIs of P
  abline(h = 0, col = "red")
  abline(v = d/(a*(e-d*h)), col = "red")
  })
legend("topright", legend = c("ZNGIs of N", "ZNGIs of P"), col = c("blue", "red"), lty = 1)

library(pracma)
peaks <- findpeaks(pop_size[, "N"])[ ,2]
periods <- peaks[length(peaks)] - peaks[length(peaks) - 1]
t_periods <- (length(times) - periods + 1):length(times)
avg_N <- mean(pop_size[t_periods, "N"])
avg_N
plot(pop_size[t_periods, "N"] ~ t_periods, type = "l")
avg_P <- mean(pop_size[t_periods, "P"])
avg_P
plot(pop_size[t_periods, "P"] ~ t_periods, type = "l")

### elliptic law
library(ggplot2)
#### The function to build Elliptic law
build_Elliptic <- function(S, C, d, sigma, rho){

  # sample coefficients in pairs
  pairs <- MASS::mvrnorm(n = S * (S-1) / 2,
                         mu = c(0, 0),
                         Sigma = sigma^2 * matrix(c(1, rho, rho, 1), 2, 2))

  # build a completely filled matrix
  M <- matrix(0, S, S)
  M[upper.tri(M)] <- pairs[,1]
  M <- t(M)
  M[upper.tri(M)] <- pairs[,2]

  # determine which connections to retain
  Connections <- (matrix(runif(S * S), S, S) <= C) * 1
  Connections[lower.tri(Connections)] <- 0
  diag(Connections) <- 0
  Connections <- Connections + t(Connections)
  M <- M * Connections

  # set diagonals
  diag(M) <- diag(M) - d
  return(M)
}

M_PP <- build_Elliptic(S = 500, C = 0.3, d = 10, sigma = 1, rho = -0.5)
M_CM <- build_Elliptic(S = 500, C = 0.3, d = 10, sigma = 1, rho = 0.5)
EVals_PP <- eigen(M_PP)$values
Re.EVals_PP <- Re(EVals_PP)
Im.EVals_PP <- Im(EVals_PP)
EVals_CM <- eigen(M_CM)$values
Re.EVals_CM <- Re(EVals_CM)
Im.EVals_CM <- Im(EVals_CM)

Re.EVals <- c(Re.EVals_PP, Re.EVals_CM)
Im.EVals <- c(Im.EVals_PP, Im.EVals_CM)
plot(Re.EVals, Im.EVals, xlab = "Real part", ylab = "Imaginary part", type = "n")
points(Re.EVals_PP, Im.EVals_PP, col = "green")
points(Re.EVals_CM, Im.EVals_CM, col = "blue")
abline(v = 0, col = "red", lty = 2)
legend("topleft", legend = c("Resource-consumer", "Competition or mutualism"), col = c("green", "blue"), pch = 1)





#### assignment A12

# niche plot
par(oma = c(2, 3, 0, 5)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(2, 2, 3, 1)) # make the plots be closer together
layout(matrix(c(1,2,3,4,5,5), ncol = 2, byrow = TRUE), heights=c(4, 4, 1))
# case 1: low niche overlap
x1 <- seq(-100, 100, by = 5)
x2 <- seq(-100, 100, by = 5) + 0.001
y1 <- dnorm(x1, mean = -50, sd = 15)
y2 <- dnorm(x1, mean = 50, sd = 15)

plot(c(-100, 100), c(0, 0.05), type = "n", ann = F, xaxt = "n" , yaxt = "n")
abline(h = 0.04, lty = 3)
points(y1 ~ x1, type = "h", col = "red")
points(y2 ~ x2, type = "h", col = "blue")
title("Low niche overlap")
text(x = -97, y = 0.05, labels = "(A)")


# case 2: high niche overlap
x1 <- seq(-50, 100, by = 5)
x2 <- seq(-50, 100, by = 5) + 0.5
y1 <- dnorm(x1, mean = 0, sd = 15)
y2 <- dnorm(x1, mean = 50, sd = 15)

plot(c(-50, 100), c(0, 0.05), type = "n", ann = F, xaxt = "n" , yaxt = "n")
abline(h = 0.04, lty = 3)
points(y1 ~ x1, type = "h", col = "red")
points(y2 ~ x2, type = "h", col = "blue")
title("High niche overlap")
text(x = -47, y = 0.05, labels = "(B)")


# case 3: different height
x1 <- seq(-50, 100, by = 5)
x2 <- seq(-50, 100, by = 5) + 0.5
y1 <- dnorm(x1, mean = 0, sd = 15) * 1.4
y2 <- dnorm(x1, mean = 50, sd = 15) * 0.5

plot(c(-50, 100), c(0, 0.05), type = "n", ann = F, xaxt = "n" , yaxt = "n")
abline(h = 0.04, lty = 3)
points(y1 ~ x1, type = "h", col = "red")
points(y2 ~ x2, type = "h", col = "blue")
text(x = -47, y = 0.05, labels = "(C)")
title("Different heights")


# case 4: a third species
x1 <- seq(-100, 100, by = 5)
x2 <- seq(-100, 100, by = 5) + 1.0
x3 <- seq(-100, 100, by = 5) + 0.5
y1 <- dnorm(x1, mean = -50, sd = 15)
y2 <- dnorm(x1, mean = 50, sd = 15)
y3 <- dnorm(x1, mean = 0, sd = 15)

plot(c(-100, 100), c(0, 0.05), type = "n", ann = F, xaxt = "n" , yaxt = "n")
abline(h = 0.04, lty = 3)
points(y1 ~ x1, type = "h", col = "red")
points(y2 ~ x2, type = "h", col = "blue")
points(y3 ~ x3, type = "h", col = "green")
text(x = -97, y = 0.05, labels = "(D)")
title("Adding a third species")


par(mai=c(0,0,0,0))
plot.new()
legend(x = "center", ncol = 2,
       legend = c("Species A", "Species B", "Species C", "Carrying capacity of resource"),
       col = c("red", "blue", "green", "black"),
       lty = c(1,1,1,3),
       bty = "n")
# axis name
mtext('Variables of resources', side = 1, outer = TRUE, line = 1)
mtext('Consumption rate/Amount of resource', side = 2, outer = TRUE, line = 1)


### figure 2
# niche plot
par(oma = c(2, 3, 0, 5)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(2, 2, 3, 1)) # make the plots be closer together
layout(matrix(c(1,1,1,2,2,2,3,3,3,3,3,3), ncol = 6, byrow = TRUE), heights=c(4, 1))
# case 1: discrete resource
x1 <- seq(-100, 100, by = 5)
x2 <- seq(-100, 100, by = 5) + 1.0
x3 <- seq(-100, 100, by = 5) + 0.5
y1 <- dnorm(x1, mean = -50, sd = 15)
y2 <- dnorm(x1, mean = 50, sd = 15)
y3 <- dnorm(x1, mean = 0, sd = 15)

plot(c(-100, 100), c(0, 0.05), type = "n", ann = F, xaxt = "n" , yaxt = "n")
abline(h = 0.04, lty = 3)
points(y1 ~ x1, type = "h", col = "red")
points(y2 ~ x2, type = "h", col = "blue")
points(y3 ~ x3, type = "h", col = "green")
title("Discrete resources")
text(x = -97, y = 0.05, labels = "(A)")


# case 2: continuous
x1 <- seq(-50, 100, by = 0.2)
x2 <- seq(-50, 100, by = 0.2) + 1.0
x3 <- seq(-50, 100, by = 0.2) + 0.5
y1 <- dnorm(x1, mean = 0, sd = 15)
y2 <- dnorm(x1, mean = 50, sd = 15)
y3 <- dnorm(x1, mean = 25, sd = 15)

plot(c(-50, 100), c(0, 0.05), type = "n", ann = F, xaxt = "n" , yaxt = "n")
abline(h = 0.04, lty = 3)
library(scales)
points(y1 ~ x1, type = "h", col = alpha("red", 0.2))
points(y2 ~ x2, type = "h", col = alpha("blue", 0.2))
points(y3 ~ x3, type = "h", col = alpha("green", 0.2))
text(x = -47, y = 0.05, labels = "(B)")
title("Continuous resources")
segments(0, 0, 0, 0.035, lty = 5)
segments(25, 0, 25, 0.035, lty = 5)
segments(40, 0, 40, 0.035, lty = 5)
arrows(0, 0.03, 25, 0.03, length = 0.05, angle = 30, code = 3)
arrows(25, 0.03, 40, 0.03, length = 0.05, angle = 30, code = 3)
text(12.5, 0.034, "d")
text(32.5, 0.034, "w")


par(mai=c(0,0,0,0))
plot.new()
legend(x = "center", ncol = 3,
       legend = c("Species A", "Species B", "Species C"),
       col = c("red", "blue", "green"),
       lty = c(1,1,1),
       bty = "n")
# axis name
mtext('Variables of resources', side = 1, outer = TRUE, line = 1)
mtext('Consumption rate/Amount of resource', side = 2, outer = TRUE, line = 1)

####
## show the there are two solution between 0 and 1
b <- seq(from = 0, to = 1, length.out = 10000)
y <- b^4 + 1 - 2*b
plot(y ~ b, type = "l")
abline(h = 0)

## find the solution between b = 0.4 to 0.6
b <- seq(from = 0.4, to = 0.6, length.out = 10000)
y <- b^4 + 1 - 2*b
plot(y ~ b, type = "l")
abline(h = 0)
ind <- which(y^2 == min(y^2))
b[ind] # 0.5437
y[ind]
