##### Part 1 #####
library(deSolve)

### Model specification
logistic_model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N*(K-N)/K # logistic growth equation
    return(list(c(dN_dt)))  # return the results
  })
}

### Model application
times <- seq(0, 10, by = 0.1)  # time steps to integrate over
state <- c(N = 10)  # initial population size
state1 <- c(N = 1000)
parms <- c(r = 1.5, K = 500)  # intrinsic growth rate and carrying capacity
parms1 <- c(r = 1.5, K = 200)

### run the ode solver
pop_size <- ode(func = logistic_model, times = times, y = state, parms = parms)
pop_size1 <- ode(func = logistic_model, times = times, y = state, parms = parms1)
pop_size2 <- ode(func = logistic_model, times = times, y = state1, parms = parms)

### plotting
plot(N ~ time, data = pop_size)
plot(N ~ time, data = pop_size1)
plot(N ~ time, data = pop_size2)

