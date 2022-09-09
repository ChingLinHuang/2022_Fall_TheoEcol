library(deSolve)
library(tidyverse)

### (1) Model specification
PP <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N - a*N*P
    dP_dt = e*a*N*P - d*P
    return(list(c(dN_dt,dP_dt)))  # return the results  
  })
}

### (2) Model application
times <- seq(0, 100, by = 0.01)  
state <- c(N = 10, P = 2)  
state1 <- c(N = 0.5/0.08 + 1, P = 1/0.1 + 1)  
parms <- c(r = 1.0, a = 0.1, e = 0.8, d = 0.5)  

# run the ode solver

pop_size <- ode(func = PP, times = times, y = state, parms = parms)
pop_size1 <- ode(func = PP, times = times, y = state1, parms = parms)

windows()
plot(P ~ N, data = pop_size1, type = "l")
# lines(P ~ N, data = pop_size, col = "red")
legend("topright", legend = c("N = 10, P = 2", "N = 10, P = 4"), 
       col = c("red", "blue"), lty = 1)

windows()
plot(P ~ time, data = pop_size, col = "red", type = "l")
lines(N ~ time, data = pop_size, col = "blue", type = "l")
legend("topright", legend = c("N", "P"), col = c("red", "blue"), lty = 1)
