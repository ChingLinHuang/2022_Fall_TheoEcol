# midterm
# 4.
# Create transition matrix
mat <- matrix(c(0   ,0   ,5.42  ,8.6   ,2.7,
                0.94,0.7 ,0.0638,0     ,0,
                0   ,0.24,0.511 ,0.0323,0,
                0   ,0   ,0.319 ,0.5484,0.093,
                0   ,0   ,0     ,0.2903,0.884), ncol = 5, byrow = T)

# Check only one dominant eigenvalue
abs(eigen(mat)$values)
# Derive dominant eigenvalue
abs(eigen(mat)$values[1])
# Derive corresponding eigenvector
as.numeric(eigen(mat)$vec[ ,1]/sum(eigen(mat)$vec[ ,1]))


initial_age <- c(100, 0, 0, 0, 0)

### for loop and matrix algebra
time <- 50
pop_size <- data.frame(Age1 = numeric(time+1),
                       Age2 = numeric(time+1),
                       Age3 = numeric(time+1),
                       Age4 = numeric(time+1),
                       Age5 = numeric(time+1))
pop_size[1, ] <- initial_age

for (i in 1:time) {
  pop_size[i+1, ] <- mat %*% as.matrix(t(pop_size[i, ]))
}

pop_size <- pop_size %>%
  round() %>%
  mutate(Total_N = rowSums(.),
         Time = 0:time) %>%
  relocate(Time)

head(round(pop_size))
pop_size[51,]/sum(pop_size[51,])


# 5.
library(deSolve)
model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dS_dt <- theta + sigma*R - beta*S*I - delta*S
    dI_dt <- beta*S*I - rho*I - gamma*I
    dR_dt <- rho*I - sigma*R - delta*R
    return(list(c(dS_dt, dI_dt, dR_dt)))  # return the results
  })
}
times <- seq(0, 10000, by = 0.01)
state <- c(S = 100, I = 2, R = 0)
parms <- c(theta = 0.1, beta = 0.01, rho = 0.05, sigma = 0, delta = 0.01, gamma = 0.02)
pop_size <- ode(func = model, times = times, y = state, parms = parms)
tail(pop_size)

# windows()
plot(S ~ time, data = pop_size[1:10000, ], type = "l", col = "red")
lines(I ~ time, data = pop_size[1:10000, ], type = "l", col = "blue")
lines(R ~ time, data = pop_size[1:10000, ], type = "l", col = "green")

# alternative
parms <- c(theta = 0.1, beta = 0.001, rho = 0.05, sigma = 0, delta = 0.01, gamma = 0.02)
pop_size <- ode(func = model, times = times, y = state, parms = parms)
tail(pop_size)
parms <- c(theta = 0.1, beta = 0.01, rho = 0.5, sigma = 0, delta = 0.01, gamma = 0.02)
pop_size <- ode(func = model, times = times, y = state, parms = parms)
tail(pop_size)

# model2 with vaccination
model2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dS_dt <- theta + sigma*R - beta*S*I - delta*S -mu*S
    dI_dt <- beta*S*I - rho*I - gamma*I
    dR_dt <- rho*I - sigma*R - delta*R
    dV_dt <- mu*S - delta*V
    return(list(c(dS_dt, dI_dt, dR_dt, dV_dt)))  # return the results
  })
}

times <- seq(0, 10000, by = 0.01)
state <- c(S = 100, I = 2, R = 0, V = 0)
parms <- c(theta = 0.1, beta = 0.01, rho = 0.05, sigma = 0, delta = 0.01, gamma = 0.02, mu = 0.005)

pop_size <- ode(func = model2, times = times, y = state, parms = parms)
tail(pop_size)
