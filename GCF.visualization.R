args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1])

savedir = args[2]
png(savedir)
hist(data[,1],
     main="Distribution fastANI results",
     sub="fraglen=30, minfrag=1, k=16",
     xlab="fraction similarity", 
     ylab="count",
     font.main=4,font.lab=7, font.sub=3,
     border="black", 
     col="blue",
     xlim=c(0,1),
     las=1, 
     breaks=60)

dev.off()
