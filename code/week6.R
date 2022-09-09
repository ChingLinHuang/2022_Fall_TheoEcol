# install.packages("deSolve")
library(deSolve)
library(tidyverse)

### (1) Model specification
PSF <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dP00_dt <- -RA*(PA0+PAA)*P00 + mA*PA0 + dA*P0A
    dPA0_dt <-  RA*(PA0+PAA)*P00 - mA*PA0 - cA*PA0
    dPAA_dt <-  cA*PA0 - mA*PAA + alpha*RA*(PA0+PAA)*P0A
    dP0A_dt <-  mA*PAA - (dA+alpha*RA*(PA0+PAA))*P0A
    return(list(c(dP00_dt, dPA0_dt, dPAA_dt, dP0A_dt)))  # return the results  
  })
}

### (2) Model application
times <- seq(0, 100, by = 0.1)  
state <- c(P00 = 0.9, PA0 = 0.09, PAA = 0.005, P0A = 0.005)  
parms <- c(RA = 0.5, mA = 0.1, cA = 0.5, dA = 0.4, alpha = 0.7)  

# run the ode solver
pop_size <- ode(func = PSF, times = times, y = state, parms = parms)

# take a look at the results
head(pop_size)
tail(pop_size)

dat <- as.data.frame(pop_size) %>% gather(key = pop, value = n, -time)
ggplot(dat, aes(x = time, y = n, color = pop)) + geom_line()
  