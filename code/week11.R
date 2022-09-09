library(tidyverse)
library(deSolve)

Prey_logistic_model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N*(1-(N/K))-a*N*P
    dP_dt = e*a*N*P-d*P
    return(list(c(dN_dt, dP_dt)))  
  })
}

times <- seq(0, 100, by = 0.01)  
state <- c(N = 40, P = 10)  
parms <- c(r = 40.0, K = 60, a = 0.1, e = 0.1, d = 0.5)  # r is chosen to be sufficiently large for time-scale separation

pop_size <- ode(func = Prey_logistic_model, times = times, y = state, parms = parms)

# population trajectories
pop_size %>%
  as.data.frame() %>%
  pivot_longer(cols = -time, names_to = "species", values_to = "N") %>%
  ggplot(aes(x = time, y = N, color = species)) + 
  geom_line(size = 1.5) +
  theme_classic(base_size = 12) +
  labs(x = "Time", y = "Population size") +
  scale_x_continuous(limits = c(0, 100.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(pop_size[, -1])*1.2), expand = c(0, 0)) +
  scale_color_brewer(name = NULL, palette = "Set1", labels = c("Prey", "Predator"), direction = -1)

# state-space diagram
pop_size %>%
  as.data.frame() %>%
  ggplot(aes(x = N, y = P)) + 
  geom_point(color = "grey60", size = 3, shape = 21) + 
  geom_vline(xintercept = with(as.list(parms), d/(e*a)), color = "#E41A1C", size = 1) +
  geom_abline(slope = with(as.list(parms), -r/(a*K)), 
              intercept = with(as.list(parms), r/a),
              color = "#377EB8", size = 1) +
  geom_point(aes(x = with(as.list(parms), d/(e*a)),
                 y = with(as.list(parms), -r/(a*K)*d/(e*a) + r/a)),
             size = 4) + 
  theme_classic(base_size = 14) +
  theme(axis.line.x = element_line(color = "#E41A1C", size = 1),
        axis.line.y = element_line(color = "#377EB8", size = 1)) +           labs(x = "Prey", y = "Predator") + 
  scale_y_continuous(limits = c(NA, 100))


# Rosenzweig-MacArthur Model
RM <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN_dt = r*N*(1-(N/K))-a*N*P/(1+a*h*N)
    dP_dt = e*a*N*P/(1+a*h*N)-d*P
    return(list(c(dN_dt, dP_dt)))  
  })
}

times <- seq(0, 2000, by = 0.01)  
state <- c(N = 40, P = 10)  
parms <- c(r = 4.0, K = 60, a = 0.1, e = 0.1, d = 0.5, h = 0.1)  


# K = 60 no preditor
pop_size <- ode(func = RM, times = times, y = state, parms = parms)
plot(pop_size[ ,2:3])

# K = 600 circle spin from outside
parms1 <- c(r = 40.0, K = 600, a = 0.1, e = 0.1, d = 0.5, h = 0.1)  
pop_size1 <- ode(func = RM, times = times, y = state, parms = parms1)
plot(pop_size1[ ,2:3])

# K = 400 circle spin from inside
state2 <- c(N = 200, P = 400)  
parms2 <- c(r = 40.0, K = 400, a = 0.1, e = 0.1, d = 0.5, h = 0.1)  
pop_size2 <- ode(func = RM, times = times, y = state2, parms = parms2)

windows()
plot(pop_size2[ ,2:3], type = "l")

# K = 250 
state3 <- c(N = 100, P = 50)  
parms3 <- c(r = 2.0, K = 250, a = 0.1, e = 0.1, d = 0.5, h = 0.1)  
pop_size3 <- ode(func = RM, times = times, y = state3, parms = parms3)

windows()
plot(pop_size3[ ,2:3], type = "l")


