# week 5

### Leslie matrix and initial age classes
leslie <- matrix(data = c(0, 1, 5,
                               0.5, 0, 0,
                               0, 0.3, 0),
                      nrow = 3,
                      ncol = 3,
                      byrow = T)

N0 <- c(10, 0, 0)

### for loop and matrix algebra
time <- 50
pop_size <- data.frame(Age1 = 0,
                       Age2 = 0,
                       Age3 = 0)
pop_size[1, ] <- N0

for (i in 2:time) {
  pop_size[i, ] <- leslie %*% t(pop_size[i-1, ])
}

# Total abundance
pop_size$N <- rowSums(pop_size)

head(pop_size)

plot(c(1,time), c(0,265), type = "n", ann = F)
lines(1:time , pop_size$Age1, col = "red")
lines(1:time , pop_size$Age2, col = "blue")
lines(1:time , pop_size$Age3, col = "green")
legend("topleft",
       legend = c("Age1", "Age2", "Age3"),
       col = c("red", "blue", "green"),
       lty = 1)


### Asymptotic growth rate and stable age distribution

asymptotic_growth <- pop_size$N[time]/pop_size$N[time-1]
asymptotic_growth

age_distribution <- pop_size[time, 1:3]/sum(pop_size[time, 1:3])
age_distribution

### Eigen-analysis of the Leslie matrix
EIGEN <- eigen(leslie)
EIGEN
as.numeric(EIGEN$values[1]) # dominant eigenvalue
as.numeric(EIGEN$vectors[, 1] / sum(EIGEN$vectors[, 1])) # corresponding eigenvector


##### Part 2
library(ggplot2)
set.seed(1234)
MAT <- matrix(rnorm(25), ncol = 5, nrow = 5)
abs(eigen(MAT)$values) # check only one dominant eigenvalue
eig_vec1 <- as.numeric(eigen(MAT)$vector[, 1])
v <- rnorm(5)
time <- 15

dat_v <- data.frame(matrix(ncol = 5, nrow = time))
dat_v[1, ] <- v
for(i in 2:time){
  dat_v[i, ] <- MAT %*% t(dat_v[i-1, ])
}

# Remake data for gganimate
dat <- data.frame(X1 = 0, X2 = 0, Time = 1)
for(i in 1:time){
  dat <- rbind(dat, data.frame(dat_v[i,1:2] / sqrt(sum(dat_v[i,1:2]^2)) * i, Time = i))
  dat <- rbind(dat, c(0,0, i+1))
}
dat <- dat[-nrow(dat), ]


ggplot(dat, aes(X1, X2, color = Time)) +
  geom_path(arrow = arrow(length = unit(0.55, "cm"))) +
  geom_abline(intercept = 0,
              slope = eig_vec1[2]/eig_vec1[1],
              color = "red",
              linetype = "dashed") # red dashed eigenvector

##### Part 3
leslie <- matrix(c(0    ,0   ,0    ,0    ,0    ,322.38,
                   0.966,0   ,0    ,0    ,0    ,0,
                   0.013,0.01,0.125,0    ,0    ,3.448,
                   0.007,0   ,0.125,0.238,0    ,30.17,
                   0.008,0   ,0    ,0.245,0.167,0.862,
                   0    ,0   ,0    ,0.023,0.75 ,0)
                 , ncol = 6, byrow = T)

eigen(leslie)$values # 2.32188

M = matrix(c(0, 0.966, 0.013, 0.007, 0.008, 0,
             0, 0, 0.01, 0, 0, 0,
             0, 0, 0.125, 0.125, 0, 0,
             0, 0, 0, 0.238, 0.245, 0.023,
             0, 0, 0, 0, 0.167, 0.75,
             322.28, 0, 3.448, 30.17, 0.862, 0), 6, 6)
max(abs(eigen(M)$value))

