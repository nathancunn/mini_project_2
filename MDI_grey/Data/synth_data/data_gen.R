set.seed(11331)
a <- mvtnorm::rmvnorm(100, mean = rep(0, 20))
b <- mvtnorm::rmvnorm(100, mean = rep(10, 20))

c <- rbind(a,b)

samp <- sample(1:nrow(c))
ind <- 1 + as.numeric(samp <= 100)

colsamp <- sample(1:ncol(c))
testdat1 <- c[ , colsamp[1:10]]
testdat2 <- c[ , colsamp[11:20]]

write.csv(ind, "MDI/Data/synth_data/ind.csv")
write.csv(testdat1, "MDI/Data/synth_data/testData1.csv")
write.csv(testdat2, "MDI/Data/synth_data/testData2.csv")


dat1 <- read.csv("MDI/Data/GaussianTestData1.csv")
hist(dat1[, 2])

dat2 <- read.csv("MDI/Data/synth_data/GaussianTestData1.csv")
hist(dat2[, 2])
