library(RColorBrewer)
library(ggplot2)
library(grid)

dataFile1 <- read.csv("MDI/Data/datagen/gaussian1.csv")
dataFile2 <- read.csv("MDI/Data/datagen/gaussian2.csv")
dataFile3 <- read.csv("MDI/Data/datagen/multinom1.csv")
dataFile4 <- read.csv("MDI/Data/datagen/multinom2.csv")


write.csv(dataFile1[, 6], "MDI/Data/datagen/gaussian_1_feature.csv")
write.csv(dataFile2[, 6], "MDI/Data/datagen/gaussian_1_feature_2.csv")


cluster <- as.numeric(dataFile1$X)

dataFile1 <- cbind(dataFile1)
gaussDF <- data.frame(rbind(dataFile1[, - 1], dataFile2[, - 1]), 
                      id = rep(c("Gaussian Test Data 1", "Gaussian Test Data 2"),
                               each = nrow(dataFile1)),
                      cluster = as.factor(rep(cluster, 2)))


ggplot(data = gaussDF) + geom_density(aes(x = f1, group = cluster, fill = cluster, alpha = 0.5)) + 
  facet_grid(.~id) + 
  scale_fill_brewer(palette = "Spectral") + 
  theme(legend.position = "none") + 
  labs(y = "Frequency", x = "")


multinomDF <- data.frame(rbind(dataFile3[, - 1], dataFile4[, - 1]), 
                      id = rep(c("Multinomial Test Data 1", "Multinomial Test Data 2"),
                               each = nrow(dataFile1)),
                      cluster = as.factor(rep(cluster, 2)))

i <- 1
out <- data.frame(cluster = rep(0, 24), id = rep(c("Multinomial Test Data 1", "Multinomial Test Data 2"),
                                                times = 12), count = rep(0, 24), val = rep(0, 24))
for(val in 1:3)
for(c in 1:4){
  for(dat in c("Multinomial Test Data 1", "Multinomial Test Data 2")) {
    out$cluster[i] <- c
    temp <- subset(multinomDF, id == dat)
    temp <- subset(temp, cluster == c)
    out$count[i] <- sum(temp$f1 == val)
    out$val[i] <- val
    i <- i + 1
  }
}
out$cluster <- as.factor(out$cluster)

ggplot(data = out, aes(group = cluster, fill = cluster, alpha = 0.5, colour = cluster)) + 
  geom_bar(aes(y = count / sum(count), x = val), stat = "identity", position = "dodge") + 
  facet_grid(.~id) + 
  scale_fill_brewer(palette = "Spectral") + 
  scale_color_manual(values = rep("black", 4))+
  theme(legend.position = "none") + 
  labs(y = "Frequency", x = "")
