rm(list=ls())
setwd("~/Bureau/SimBioProj/r_analysis/")
#install.packages("varhandle")
library(varhandle)
library(ggplot2)

# read file which contain the historic of fintess and mutations
f = read.table('fitness_by_event_500.txt', sep =' ', header=F, dec='.')
ev <- c()
fit <- c()
for (i in 1:length(f[2,])){
  fit[i] <- as.numeric(unfactor(f[2,i]))
  ev[i] <- as.character(f[1,i])
}

# fitness change
ev <- ev[-1]
fit <- fit[-length(fit)] - fit[-1]
par(mfrow=c(1,1))
hist(fit, col="red3", nclass=20)

# using only mutations kept (changing fitness)
data <- cbind.data.frame(ev,fit)
names(data) <- c("e","f")
data <- data[(data$f > 1e-10 | data$f < -1e-10), ]

# see the change of fitness induced by each mutation
par(mfrow=c(2,2))
hist(data$f, main='all', xlab='Fitness Change', sub=paste('mean =', signif(mean(data$f),4)), col="red3", nclass=20)
hist(data$f[data$e=="del"], main='del', xlab='Fitness Change', sub=paste('mean =', signif(mean(data$f[data$e=="del"]),4)), col="red3", nclass=20)
hist(data$f[data$e=="ins"], main='ins', xlab='Fitness Change', sub=paste('mean =', signif(mean(data$f[data$e=="ins"]),4)), col="red3", nclass=20)
hist(data$f[data$e=="inv"], main='inv', xlab='Fitness Change', sub=paste('mean =', signif(mean(data$f[data$e=="inv"]),4)), col="red3", nclass=20)

# global view to compare the mutations
ggplot(data,aes(x=data$f,group=data$e,fill=data$e))+
  geom_histogram(position="dodge",binwidth=5e-05)+theme_bw()

