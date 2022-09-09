# week 4

library(tidyverse)

### (1) Set the parameters
r <- 2.9
K <- 500
N0 <- 11
time <- 100

### (2) Define the discrete logistic growth equation
log_fun <- function(r, N, K){N + r*N*(1-N/K)}  

### (3) Use for loop to iterate over the time sequence
pop_size <- numeric(time)
pop_size[1] <- N0

for (i in 2:time) {pop_size[i] <- log_fun(r = r, N = pop_size[i - 1], K = K)}

pop_data <- pop_size %>% 
  as.data.frame() %>% 
  rename(., pop_size = `.`) %>%
  mutate(time = 0:(time-1)) %>%
  relocate(time)

head(pop_data)

##   time  pop_size
## 1    0  10.00000
## 2    1  27.64000
## 3    2  74.64171
## 4    3 188.93980
## 5    4 400.51775
## 6    5 543.95762

  
  ### Population trajectory
  ggplot(pop_data, aes(x = time, y = pop_size)) + 
  geom_point() + 
  geom_line() +
  geom_hline(yintercept = K, color = "red", size = 1.2, linetype = "dashed") + 
  geom_text(x = time*1.02, y = K+50, label = "italic(K)", color = "red", size = 6.5, parse = T) +
  labs(y = expression(italic(N)), title = paste0("Discrete logistic growth", "\n", "(r = ", r, ", K = ", K, ", N0 = ", N0, ")")) + 
  scale_x_continuous(limits = c(0, time*1.05), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, max(pop_size)*1.1), expand = c(0, 0)) + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))

### Cobweb plot/logistic map
cobweb_data <- data.frame(Nt = rep(pop_size[-time], each = 2)[-1], 
                          Nt1 = rep(pop_size[-1], each = 2)[-length(rep(pop_size[-1], each = 2))])

#  Nt   Nt1
#----------
#  N1   N2
#  N2   N2
#  N2   N3

logistic_map <- data.frame(Nt = seq(0, (r+1)/r*K, by = 0.1)) %>%
  mutate(Nt1 = Nt + r*Nt*(1-Nt/K))

ggplot() + 
  geom_line(data = logistic_map, aes(x = Nt, y = Nt1), color = "green", size = 1.2) + 
  geom_path(data = cobweb_data, aes(x = Nt, y = Nt1), color = "blue", size = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1) + 
  labs(x = expression(italic(N[t])),
       y = expression(italic(N[t+1])), 
       title = paste0("Cobweb plot/logistic map", "\n", "(r = ", r, ", K = ", K, ", N0 = ", N0, ")")) + 
  scale_x_continuous(limits = c(0, (r+1)/r*K*1.05), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, max(pop_size)*1.1), expand = c(0, 0)) + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
