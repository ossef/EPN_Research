

modelname1 = "Courbes.res"
data1 = read.table(modelname1)
attach(data1);


 # library
library(ggplot2)
 

    
ggplot(data1, aes(V2))+
  geom_histogram()
ggplot(data1, aes(V3))+
  geom_histogram()


#barplot(H, names.arg=P ,xlab = "Number of panels", ylim=c(0,6), ylab="Mean waiting time of a DP",border="red",col="blue", density=14,beside=T)
#barplot(Time, names.arg=Experience ,xlab = "Experiment", ylim=c(0,1.5), ylab="Mean waiting time of a DP in the EPN",border="red",col="blue", density=14)
#barplot(Loss, names.arg=Experience ,xlab = "Experiment", ylim=c(0,35), ylab="Loss rate of EPs in the EPN", border="red",col="blue", density=14)
