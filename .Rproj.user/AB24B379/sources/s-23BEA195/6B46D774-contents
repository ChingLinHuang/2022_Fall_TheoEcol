library(deSolve)
library(tidyverse)

### Model specification
PSF <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dP00_dt = P0A*dA + PA0*mA - P00*(PA0 + PAA)*rA
    dPA0_dt = P00*(PA0 + PAA)*rA - PA0*mA - PA0*cA
    dPAA_dt = PA0*cA - PAA*mA + P0A*(PA0 + PAA)*rA*alpha
    dP0A_dt = PAA*mA - P0A*(PA0 + PAA)*rA*alpha - P0A*dA

    return(list(c(dP00_dt, dPA0_dt, dPAA_dt, dP0A_dt)))
  })
}

### Model parameters
times <- seq(0, 20, by = 0.1)
state <- c(P00 = 0.25, PA0 = 0.25, PAA = 0.25, P0A = 0.25)
parms <- c(rA = 0.5, mA = 0.1, cA = 0.5, dA = 0.4, alpha = 0.7)

### ODE solver
pop_size <- ode(func = PSF, times = times, y = state, parms = parms)

# take a look at the results
head(pop_size)
tail(pop_size)

### Visualization
pop_size %>%
  as.data.frame() %>%
  gather(key = "patch", value = "proportion", -time) %>%
  ggplot(aes(x = time, y = proportion, color = patch)) +
  geom_line(size = 1.5)

plot(range(times), c(0,1), type = "n", xlab = "time", ylab = "proportion")
lines(P00 ~ time, data = pop_size, col = "tomato")
lines(P0A ~ time, data = pop_size, col = "navy")
lines(PA0 ~ time, data = pop_size, col = "gray")
lines(PAA ~ time, data = pop_size, col = "orange")
legend("topleft", legend = c("P00", "P0A", "PA0", "PAA"), col = c("tomato", "navy", "gray", "orange"), lty = 1)


#####  Assignment 6
library(deSolve)

MP_model <- function(c_0, e_0){

  MP <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dP_dt = c_0*P*(1 - h - P) - e_0*P

      return(list(dP_dt))
    })
  }

  # dataframe for saving final fate of population size with corresponding h
  h_seq <- seq(0, 1, by = 0.01)
  dat <- data.frame(h = h_seq, P = 0)
  for (ind in 1:length(h_seq)){
    times <- seq(0, 500, by = 0.1)
    state <- c(P = 0.5)
    parms <- c(c_0 = c_0, e_0 = e_0, h = h_seq[ind])

    pop_size <- ode(MP, times = times, y = state, parms = parms)

    dat$P[ind] <- pop_size[length(times), 2]
  }
  plot(P ~ h, data = dat, type = "l", main = paste0("e = ", e_0, ", c = ", c_0))
  abline(v = 1 - e_0/c_0, lty = 2)
}

MP_model(c_0 = 0.1, e_0 = 0.02)
MP_model(c_0 = 0.1, e_0 = 0.05)
MP_model(c_0 = 0.1, e_0 = 0.07)




