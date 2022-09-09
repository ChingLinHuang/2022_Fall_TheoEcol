# week 5

library(tidyverse)

#### Part 1
### Leslie matrix and initial age classes
leslie_mtrx <- matrix(data = c(0, 1, 5,
                               0.5, 0, 0,
                               0, 0.3, 0),
                      nrow = 3, 
                      ncol = 3,
                      byrow = T)

initial_age <- c(10, 0, 0)

### for loop and matrix algebra
time <- 50
pop_size <- data.frame(Age1 = numeric(time+1),
                       Age2 = numeric(time+1),
                       Age3 = numeric(time+1))
pop_size[1, ] <- initial_age

for (i in 1:time) {
  pop_size[i+1, ] <- leslie_mtrx %*% as.matrix(t(pop_size[i, ]))
}

pop_size <- pop_size %>% 
  round() %>%
  mutate(Total_N = rowSums(.), 
         Time = 0:time) %>%
  relocate(Time)

head(round(pop_size)) 

### Asymptotic growth rate and stable age distribution 
asymptotic_growth <- round(pop_size[time+1, 5]/pop_size[time, 5], 3)
asymptotic_growth

age_distribution <- round(pop_size[time+1, 2:4]/sum(pop_size[time+1, 2:4]), 3)
age_distribution

### Eigenanalysis of the Leslie matrix
eigen_out <- eigen(leslie_mtrx)
as.numeric(eigen_out$values[1]) %>% round(., 3)  # dominant eigenvalue
as.numeric(eigen_out$vectors[, 1]/sum(eigen_out$vectors[, 1])) %>% 
  round(., 3)  # stable age distribution

#### Part 2 - Visualizing population dynamics and age distribution
### Population sizes for each age class
pop_size %>%
  pivot_longer(cols = -Time, names_to = "Age_class", values_to = "N") %>%
  ggplot(aes(x = Time, y = N, color = Age_class)) + 
  geom_point() + 
  geom_line() + 
  labs(x = "time", y = expression(italic(N))) +
  theme_classic(base_size = 12) +
  scale_x_continuous(limits = c(0, time*1.05), expand = c(0, 0)) +
  scale_y_log10(limits = c(1, max(pop_size$Total_N)*1.05), expand = c(0, 0)) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "black"),
                     name = NULL,
                     label = c("Age1", "Age2", "Age3", "Total"))

### Stable age distribution
library(gganimate)

age_animate <- pop_size %>% 
  mutate(across(Age1:Age3, function(x){x/Total_N})) %>%
  select(Time, Age1:Age3) %>%
  pivot_longer(Age1:Age3, names_to = "Age", values_to = "Proportion") %>%
  ggplot(aes(x = Age, y = Proportion, fill = Age)) + 
  geom_bar(stat = "identity", show.legend = F) +
  labs(x = "") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1") + 
  theme_classic(base_size = 12) + 
  transition_manual(Time) + 
  ggtitle("Time {frame}") + 
  theme(title = element_text(size = 15))

anim_save("age_distribution.gif", age_animate, nframes = time + 1, fps = 4, width = 5, height = 4, units = "in", res = 300)


##### Part 3
leslie <- matrix(c(0,0,0,0,0,322.38,0.966,0,0,0,0,0,0.013,0.01,0.125,0,0,3.448,0.007,0,0.125,0.238,0,30.17,0.008,0,0.023,0.245,0.167,0.862,0,0,0,0,0.75,0), ncol = 6, byrow = T)

eigen(leslie)$values[1]

