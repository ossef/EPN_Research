modelname1 = "Hybrid.res"
data1 = read.table(modelname1)
attach(data1);

Experience = V1
Lambda1 = V2
Lambda2 = V3
Lambda3 = V4
Lambda4 = V5
Lambda5 = V6
Lambda6 = V7
Lambda7 = V8
Number = V9
Time = V10
Loss = V11


barplot(Number, names.arg=Experience ,xlab = "Experiment", ylim=c(0,45), ylab="Mean Number of DPs in the EPN",border="red",col="blue", density=14)
barplot(Time, names.arg=Experience ,xlab = "Experiment", ylim=c(0,1.5), ylab="Mean waiting time of a DP in the EPN",border="red",col="blue", density=14)
barplot(Loss, names.arg=Experience ,xlab = "Experiment", ylim=c(0,35), ylab="Loss rate of EPs in the EPN", border="red",col="blue", density=14)
