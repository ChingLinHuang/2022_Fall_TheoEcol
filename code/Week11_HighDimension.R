########################################################################################################################
########################################################################################################################
#### Teaching material: Theoretical Ecology ############################################################################
#### High-dimension systems and random matrix ##########################################################################
########################################################################################################################
########################################################################################################################



########################################################################################################################
########################################################################################################################
#### Set your working directory and load data files
########################################################################################################################
########################################################################################################################
setwd('/Users/Chris/Dropbox/Public/TEACHING/IntroductionTheoreticalEcology')



########################################################################################################################
########################################################################################################################
#### Install the required R packages
########################################################################################################################
########################################################################################################################
#### load the package into active memory
require('deSolve')
require('ggplot2')
require('tidyr')
require('dplyr')



# ########################################################################################################################
# ########################################################################################################################
# #### Generalized Lotka-Volterra model ##################################################################################
# ########################################################################################################################
# ########################################################################################################################
# #### Lotka-Volterra model
# GLV = function(Time, State, Pars){
#   with(as.list(c(State, Pars)),{
#     # dNX = State * (r + A %*% State)
#     dNX = diag(State) %*% (r + A %*% State)
#     return(list(c(dNX)))
#   })
# }
#
# #### Parameter setting 1 -- Fig1a
# r = rep(1, 3)
# A = -matrix(c(10, 9, 5,
#               9, 10, 9,
#               5, 9, 10), 3, 3, byrow = TRUE)
# parms = list(r = r, A = A)
# times = seq(0, 50, by = 0.1)
# ini= c(0.2, 0.1, 0.1)
#
# #### Parameter setting 2 -- Fig1b
# r = rep(1, 3)
# A = -matrix(c(10, 7, 12,
#               15, 10, 8,
#               7, 11, 10), 3, 3, byrow = TRUE)
# parms = list(r = r, A = A)
# times = seq(0, 800, by = 0.1)
# ini= c(0.2, 0.1, 0.1)
#
# #### Parameter setting 3 -- Fig2a
# r = rep(1, 3)
# A = -matrix(c(10, 6, 12,
#               14, 10, 2,
#               8, 18, 10), 3, 3, byrow = TRUE)
# parms = list(r = r, A = A)
# times = seq(0, 200, by = 0.1)
# ini= c(0.05, 0.05, 0.05)
#
# #### Parameter setting 3 -- Fig2b
# r = c(1, 0.72, 1.53, 1.27)
# A = -matrix(c(1, 1.1, 1.5, 0,
#               0, 0.7, 0.3, 0.9,
#               3.5, 0, 1.5, 0.7,
#               1.5, 0.6, 0.4, 1.3), 4, 4, byrow = TRUE)
# parms = list(r = r, A = A)
# times = seq(0, 200, by = 0.1)
# ini= c(0.05, 0.05, 0.05, 0.05)
#
# #### Run the simulation
# out = as.data.frame(ode(func=GLV, y=ini, parms=parms, times=times))
#
# #### Gather data for ploting
# comp.gg = gather(out, key=state.var, value=count, -time)
#
# #### Plot
# ggplot(comp.gg,
#        aes(x=time, y=count, color=state.var)) +
#   geom_line() +
#   geom_hline(yintercept=0, linetype=1, color="black") +
#   # scale_x_continuous(expand=c(0, 0)) +
#   scale_y_continuous(expand=c(0, 0)) +
#   ylab(expression(paste("Population size"))) +
#   xlab("Time") +
#   theme_bw() +
#   theme(legend.position="none",
#         plot.title=element_text(hjust=0.5),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.background=element_rect(colour="black"))
#
# #### Solve analytical solution
# solve(A, -r)
# eigen(A)$values
#
# #### Get average abundance for cycles and chaos
# Data =
#   comp.gg %>% filter(time > 50) %>%
#   group_by(state.var) %>%
#   summarize(average = mean(count))
# Data$time = 200
#
# #### Plot with average and analytical abundance
# ggplot(comp.gg,
#        aes(x=time, y=count, color=state.var)) +
#   geom_line() +
#   geom_hline(yintercept=0, linetype=1, color="black") +
#   geom_hline(data=Data, aes(yintercept=average, color=state.var), linetype=3) +
#   geom_point(data=Data, aes(x=time, y=average, color=state.var), size=2) +
#   # scale_x_continuous(expand=c(0, 0)) +
#   scale_y_continuous(expand=c(0, 0)) +
#   ylab(expression(paste("Population size"))) +
#   xlab("Time") +
#   theme_bw() +
#   theme(legend.position="none",
#         plot.title=element_text(hjust=0.5),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.background=element_rect(colour="black"))





########################################################################################################################
########################################################################################################################
#### May's random matrix approach ######################################################################################
########################################################################################################################
########################################################################################################################
#### The function to build May's random matrix
build_May <- function(S, C, d, mean = 0, sigma){

  # fill the whole matrix
  M <- matrix(rnorm(S * S, mean = mean, sd = sigma), nrow = S, ncol = S)

  # remove connections
  remove <- matrix(runif(S * S) <= C, nrow = S, ncol = S)
  M <- M * remove

  # set diagonals
  diag(M) <- diag(M) - d

  return(M)
}


#### Function to plot eigenvalues
plot_eigenvalues <- function(M, prediction = NULL){
  eig <- eigen(M)$values
  data <- data.frame(Real = Re(eig),
                     Imaginary = Im(eig))
  plot(Imaginary ~ Real, data)
  abline(v = 0, col="red", lty=3)
  abline(v = prediction)
}


#### Run
N <- 500
connectance <- 0.3
self.limit <- 10
variance <- 1
M <- build_May(S=N, C=connectance, d=self.limit, sigma=variance)
plot_eigenvalues(M)


#### The function to build Elliptic law
build_Elliptic <- function(S, C, d, sigma, rho){

  # sample coefficients in pairs
  pairs <- MASS::mvrnorm(n = S * (S-1) / 2,
                         mu = c(0, 0),
                         Sigma = sigma^2 * matrix(c(1, rho, rho, 1), 2, 2))

  # build a completely filled matrix
  M <- matrix(0, S, S)
  M[upper.tri(M)] <- pairs[,1]
  M <- t(M)
  M[upper.tri(M)] <- pairs[,2]

  # determine which connections to retain
  Connections <- (matrix(runif(S * S), S, S) <= C) * 1
  Connections[lower.tri(Connections)] <- 0
  diag(Connections) <- 0
  Connections <- Connections + t(Connections)
  M <- M * Connections

  # set diagonals
  diag(M) <- diag(M) - d
  return(M)
}


#### Function to plot eigenvalues ggplot version
plot_eigenvalues <- function(M){
  eig <- eigen(M)$values
  data <- data.frame(Real = Re(eig),
                     Imaginary = Im(eig))
  p <- ggplot(data) + aes(x = Real, y = Imaginary) +
    geom_point() +
    coord_equal() +
    geom_vline(xintercept = 0, colour = "red", linetype = 2)
  show(p)
}


#### Run Elliptic law
N <- 500
connectance <- 0.3
self.limit <- 10
variance <- 1
correlation <- -0.5
M <- build_Elliptic(S=N, C=connectance, d=self.limit, sigma=variance, rho=correlation)
plot_eigenvalues(M)


# #### Generate random matrix with off-diagonal elements having non-zero means
# N = 500
# connectance = 0.7
# self.limit = 10
# variance = 1
# average = -0.1
# M = build_May(n=N, C=connectance, d=self.limit, mean=average, sigma=variance)
# plot_eigenvalues(M)
# abline(v = (-self.limit) + (N - 1) * connectance * average)


# #### Generate modular random matrix
# N = 500
# connectance1 = 0.7
# connectance2 = 0.5
# self.limit = 10
# variance = 1
# average1 = 0.1
# average2 = 0.2
# M1 = build_May(n=N, C=connectance1, d=self.limit, mean=average1, sigma=variance)
# M2 = build_May(n=N, C=connectance2, d=self.limit, mean=average2, sigma=variance)
# M3 = build_May(n=N, C=0.05, d=0, mean=0, sigma=variance)
# M4 = rbind(cbind(M1, M3), cbind(M3, M2))
# plot_eigenvalues(M4)
# abline(v = (-self.limit) + (N - 1) * connectance1 * average1)
# abline(v = (-self.limit) + (N - 1) * connectance2 * average2)
