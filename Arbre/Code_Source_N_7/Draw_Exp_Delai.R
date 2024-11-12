
modelname1 = "Test.res"
data1 = read.table(modelname1)
attach(data1);

Experience = V1
delay = V11



barplot(delay, names.arg=Experience ,ylim =c(0,5),xlab = "Different configurations with 30 panels", ylab="Mean waiting time of a DP in the EPN",border="red",col="blue", density=14)
