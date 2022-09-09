###### part 1 ######
# install.packages("deSolve")
library(deSolve)

### (1) Model specification
exponential_model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N  # exponential growth equation
    return(list(c(dN_dt)))  # return the results
  })
}

### (2) Model application
times <- seq(0, 10, by = 0.1)  # time steps to integrate over
state <- c(N = 10)  # initial population size
parms <- c(r = 1.5)  # intrinsic growth rate

# run the ode solver
pop_size <- ode(func = exponential_model, times = times, y = state, parms = parms)

# take a look at the results
head(pop_size)

### (3) Plotting
plot(N ~ time, data = pop_size)
plot(N ~ time, data = pop_size, log = "y")

###### part 2 ######

times <- seq(0, 10, by = 0.1)  # time steps to integrate over
state <- c(N = 10)  # initial population size
parms <- c(r = 1.5)  # intrinsic growth rate
# default: LSODA
pop_size <- ode(func = exponential_model, times = times, y = state, parms = parms)

# Euler's method
pop_size_1 <- ode(func = exponential_model, times = times, y = state, parms = parms, method = "euler")

# compare different method
par(mfrow = c(1,2))
plot(N ~ time, data = pop_size, main = "LSODA")
curve(state[1]*exp(parms[1]*x), times[1], times[length(times)], col = "red", add = T) # correct curve
plot(N ~ time, data = pop_size_1, main = "Euler")
curve(state[1]*exp(parms[1]*x), times[1], times[length(times)], col = "red", add = T) # correct curve

# compare different method in log transformation
par(mfrow = c(1,2))
plot(N ~ time, data = pop_size, log = "y", main = "LSODA")
curve(state[1]*exp(parms[1]*x), times[1], times[length(times)], col = "red", add = T) # correct curve
plot(N ~ time, data = pop_size_1, log = "y", main = "Euler")
curve(state[1]*exp(parms[1]*x), times[1], times[length(times)], col = "red", add = T) # correct curve


###### part 3 ######
### Model specification
exponential_model_fluc <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = (r_hat + sigma*sin(omega*times + phi)/2)*N  # exponential growth equation
    return(list(c(dN_dt)))  # return the results
  })
}
### Parameters
times <- seq(0, 10, by = 0.1)  # time steps to integrate over
state <- c(N = 10)  # initial population size
parms <- c(r_hat = 1.5, sigma = 1, omega = 2*pi, phi = 0)  # intrinsic growth rate

### Fluctuating growth rate
r = parms[1] + parms[2]/2*sin(parms[3]*times + parms[4])
plot(r ~ times, type = "l")

### Solving model
pop_size <- ode(func = exponential_model_fluc, times = times, y = state, parms = parms)

### Plotting
plot(N ~ times, data = pop_size)
curve(state[1]*exp(parms[1]*x - parms[2]/2/parms[3]*(cos(parms[3]*x + parms[4]) - cos(parms[4]))), add = T) # correct curve
plot(N ~ times, data = pop_size, log = "y")
curve(state[1]*exp(parms[1]*x - parms[2]/2/parms[3]*(cos(parms[3]*x + parms[4]) - cos(parms[4]))), add = T) # correct curve
