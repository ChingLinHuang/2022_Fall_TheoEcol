library(ggplot2)
library(tidyverse)
library(deSolve)
library(gganimate)
library(gifski)

time <- seq(0, 10, by = 0.001)
e1 <-


  ### Model specification
  ERROR <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      de1_dt = c1 * e1 + c2 * e2
      de2_dt = d1 * e1 + d2 * e2
      return(list(c(de1_dt, de2_dt)))
    })
  }

  ### Model parameters
  times <- seq(0, 5, by = 0.0001)
  state <- c(e1 = 1, e2 = 1)
  parms <- c(c1 = -1, c2 = 1, d1 = -2, d2 = -1)

  ### Model application
  error <- ode(func = ERROR, times = times, y = state, parms = parms)
  tail(error)

  plot(e1 ~ e2, error, type = "l")

  ### Visualize the population dynamics
p <- error %>%
    as.data.frame() %>%
    ggplot(aes(x = e1, y = e2)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype="dashed", color = "red") +
    geom_hline(yintercept = 0, linetype="dashed", color = "red") +
    labs(subtitle = "Time: {round(frame_time, digit = 1)}") +
    transition_time(time) +
    shadow_wake(wake_length = 1)
gif <- animate(p, renderer = gifski_renderer())
setwd("C:\\Users\\andyh\\OneDrive\\Documents\\2022_Fall_TheoEcol\\assignment_figures")
anim_save(filename = "W9_dynamics_error.gif", gif)
