library(deSolve)

# 2.(7)
PD <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dP1_dt = c1*P1*(1-P1) - m*P1
    dP2_dt = c2*P2*(1-P1-P2) - m*P2 - c1*P1*P2
    dP3_dt = c3*P3*(1-P1-P2-P3) - m*P3 - c1*P1*P3 - c2*P2*P3
    return(list(c(dP1_dt, dP2_dt, dP3_dt)))  
  })
}


times <- seq(0, 2000, by = 0.01)  
state <- c(P1 = 0.1, P2 = 0.1, P3 = 0.1)  
parms <- c(c1 = 0.2, c2 = 0.5, c3 = 0.7, m = 0.1)

pop_size <- ode(func = PD, times = times, y = state, parms = parms)
tail(pop_size)

# 4.(2)-
Succ <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    df_dt = m*(S+E+M+R) - cl*(S+M+R)*f - ce*(E+M)*f
    dE_dt = ce*(E+M)*f - m*E - cl*(S+M+R)*E
    dS_dt = cl*(S+M+R)*f - m*S- ce*(E+M)*S - r*S
    dM_dt = ce*(E+M)*S + cl*(S+M+R)*E - m*M - r*M
    dR_dt = r*(S+M) - m*R
    return(list(c(df_dt, dE_dt, dS_dt, dM_dt, dR_dt)))  
  })
}

times <- seq(0, 2000, by = 0.01)  
state <- c(f = 0.96, E = 0.03, S = 0.01, M = 0, R = 0)  
parms <- c(ce = 0.8, cl = 0.1, r = 0.8, m = 0.05)

#(2)
pop_size <- ode(func = Succ, times = times, y = state, parms = parms)
tail(pop_size)

#(3)
parms1 <- c(ce = 0.8, cl = 0.1, r = 0.08, m = 0.05)

pop_size1 <- ode(func = Succ, times = times, y = state, parms = parms1)
tail(pop_size1)

#(4)
parms2 <- c(ce = 0.8, cl = 0.1, r = 0.8, m = 0.075)

pop_size2 <- ode(func = Succ, times = times, y = state, parms = parms2)
tail(pop_size2)