##### Part 1: Model the discrete logistic population growth using for loops
### (1) Define the discrete logistic growth equation
log_fun <- function(r, N, K){N + r*N*(1-N/K)}

### (2) Set the parameters
r <- 2
K <- 500
N0 <- 10
time <- 100

### (3) Use for loop to iterate over the time sequence
pop_size <- data.frame(times = 1:time)
pop_size$N[1] <- N0
head(pop_size)

for(i in 2:time){
  pop_size$N[i] <- log_fun(r = r, N = pop_size$N[i - 1], K = K)
}

head(pop_size)

### (4) Population trajectory
plot(N ~ times, data = pop_size, type = "l")
abline(h = K, col = "red")
points(N ~ times, data = pop_size)


###### Part 2: Bifurcation curve
### (1) data setting
# intrinsic growth rate sequence
r_seq <- seq(from = 2.82, to = 2.83, by = 0.001)
# number of sampling
N_rep <- 5000

# data
dat_plot <- data.frame(r = rep(r_seq, each = N_rep), N = 0)
head(dat_plot)

### (2) Run the discrete logistic model
for (lo in 1:length(r_seq)){
  log_fun <- function(r, N, K){N + r*N*(1-N/K)}

  r <- r_seq[lo]
  K <- 500
  N0 <- 10
  time <- 50000

  pop_size <- data.frame(times = 1:time)
  pop_size$N[1] <- N0

  for(i in 2:time){
    pop_size$N[i] <- log_fun(r = r, N = pop_size$N[i - 1], K = K)
  }

  # save the data
  dat_plot$N[(1 + (lo - 1)*N_rep):(lo*N_rep)] <- pop_size$N[(nrow(pop_size) - N_rep + 1):nrow(pop_size)]
}

plot(N ~ r, data = dat_plot, cex = 0.1, pch = 20)
