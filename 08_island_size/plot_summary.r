x <- read.table("summary.txt", sep="\t", stringsAsFactors=F, header=T, comment="")

# modify scale of some numbers
#ROH in Mbp
x[,7] <- x[,7] / 1000000
# het x 100
x[,8] <- x[,8] * 100

# plot all
par(mfrow=c(2,3))
par(mar=c(2, 4.5, 1, 1))
plot(x[,5], x[,12], pch=19, xlab="", xaxt="n", ylab="MSMC Pop. Size (20 Kya)", cex.lab=1.1, cex.axis=1.1)
axis(1, at=c(4,5,6,7,8), labels=F)
abline(lm(x[,12] ~ x[,5]))
summary(lm(x[,12] ~ x[,5]))

plot(x[,5], x[,8], pch=19, xlab="", xaxt="n", ylab="Heterozygosity (x 100)", cex.lab=1.1, cex.axis=1.1)
axis(1, at=c(4,5,6,7,8), labels=F)
abline(lm(x[,8] ~ x[,5]))
summary(lm(x[,8] ~ x[,5]))

plot(x[,5], x[,9], pch=19, xlab="", xaxt="n", ylab="Std. Dev. Heterozygosity", cex.lab=1.1, cex.axis=1.1)
axis(1, at=c(4,5,6,7,8), labels=F)
abline(lm(x[,9] ~ x[,5]))
summary(lm(x[,9] ~ x[,5]))

plot(x[,5], x[,11], pch=19, xlab="Island Size (log)", ylab="# Homozygous Non-Reference ERVs", cex.lab=1.1, cex.axis=1.1)
abline(lm(x[,11] ~ x[,5]))
summary(lm(x[,11] ~ x[,5]))

plot(x[,5], x[,7], pch=19, xlab="Island Size (log)", ylab="Length ROH (Mbp)", cex.lab=1.1, cex.axis=1.1)
abline(lm(x[,7] ~ x[,5]))
summary(lm(x[,7] ~ x[,5]))

plot(x[,5], x[,6], pch=19, xlab="Island Size (log)", ylab="# ROH", cex.lab=1.1, cex.axis=1.1)
abline(lm(x[,6] ~ x[,5]))
summary(lm(x[,6] ~ x[,5]))

# without luteirostris
x2 <- x[-3,]

summary(lm(x2[,7] ~ x2[,5]))
summary(lm(x2[,6] ~ x2[,5]))

# plot separate two plots with just two different ways of estimating pop size
par(mfrow=c(1,2))
par(mar=c(4.5, 4.5, 1, 1))
plot(x[,5], x[,12], pch=19, xlab="Island Size (log)", ylab="MSMC Pop. Size (20 Kya)", cex.lab=1.1, cex.axis=1.1)
abline(lm(x[,12] ~ x[,5]))
summary(lm(x[,12] ~ x[,5]))
plot(x[,5], x[,13], pch=19, xlab="Island Size (log)", ylab="MSMC Harmonic Mean Pop. Size", cex.lab=1.1, cex.axis=1.1)
abline(lm(x[,13] ~ x[,5]))
summary(lm(x[,13] ~ x[,5]))




