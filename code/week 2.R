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
# original setting
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

# minimize the time step
times <- seq(0, 10, by = 0.01)  # time steps to integrate over
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


###### part 3 ######
### Model specification
exponential_model_fluc <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = (r_hat + sigma*sin(omega*times))*N  # exponential growth equation
    return(list(c(dN_dt)))  # return the results
  })
}
### Parameters
times <- seq(0, 10, by = 0.1)  # time steps to integrate over
state <- c(N = 10)  # initial population size
parms <- c(r_hat = -0.1, sigma = 5, omega = 2*pi)  # intrinsic growth rate

### Fluctuating growth rate
r = parms[1] + parms[2]*sin(parms[3]*times)
plot(r ~ times, type = "l")

### Solving model
pop_size <- ode(func = exponential_model_fluc, times = times, y = state, parms = parms)

### Plotting
plot(N ~ times, data = pop_size)
curve(state[1]*exp(parms[1]*x - parms[2]/parms[3]*(cos(parms[3]*x) - 1)), add = T, col = "red") # correct curve
plot(N ~ times, data = pop_size, log = "y")
curve(state[1]*exp(parms[1]*x - parms[2]/parms[3]*(cos(parms[3]*x) - 1)), add = T, col = "red") # correct curve



########### A2
library(deSolve)

exponential_model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N + I
    return(list(c(dN_dt)))})
}
times <- seq(0, 7, by = 0.1)
state <- c(N = 10)
parms <- c(r = 1.2, I = 3) # parameters
parms_no_I <- c(r = 1.2, I = 0)

pop_size <- ode(func = exponential_model, times = times,
                y = state, parms = parms)
pop_size_no_I <- ode(func = exponential_model, times = times,
                     y = state, parms = parms_no_I)

# plot population dynamics
plot(N ~ time, data = pop_size, type = "l", col = "red")
lines(N ~ time, data = pop_size_no_I, type = "l", col = "blue")
legend("topleft",
       legend = c("with immigration", "without immigration"),
       col = c("red", "blue"),
       lty = 1)

dif <- pop_size[ ,"N"] - pop_size_no_I[ ,"N"]
plot(times, dif, type = "l")
