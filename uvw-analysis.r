setwd("I:\\Workspace\\windsim\\windspeed")

library(sandwich)
library(forecast)
library(fitdistrplus)
library(tseries)
xmat <- read.table("uvwp6.txt", sep = ",", skip = 4)

xmat.orig <- xmat


###10 minute nobs
nmin <- 2
nobs <- 100*60*nmin

 start.obs <- 195000
end.obs <- start.obs + nobs

start.obs <- 170000
end.obs  <- 200000

xmat <- xmat.orig[start.obs:end.obs,]

u.low <- xmat[,2]
v.low <- xmat[,3]
w.low <- xmat[,4]

u.mid <- xmat[,6]
v.mid <- xmat[,7]
w.mid <- xmat[,8]

u.hi <- xmat[,10]
v.hi <- xmat[,11]
w.hi <- xmat[,12]

mag.hi <- sqrt(u.hi**2 + v.hi**2)

mag.spline <- smooth.spline(mag.hi, df = 10)
plot(mag.hi, type = "l")
lines(mag.spline, col = "red")
mag.resids <- mag.hi - mag.spline$y

mag.sub <- mag.resids[10000:14000]
plot(mag.sub, type = "l")
adf.test(mag.sub)
manual.arima <- Arima(mag.resids, order = c(2, 1, 1))
acf(manual.arima$residuals)
mag.sim.man <- arima.sim(manual.arima, 10000)

mag.arima <- auto.arima(mag.sub)
mag.resids <- mag.arima$residuals
acf(mag.resids)
mag.arima$coef
arvec <- mag.arima$arma
minroots <- min(Mod(polyroot(c(1, - coef(mag.arima)))))

asim <- arima.sim(n = 5000, list(ar = c(0.2929, -0.1570), 
          sd = sqrt(0.1493)))

plot(asim, type = "l")

ts.sim <- arima.sim(list(order = c(2,1,0), ar =c(0.2929  ,  -0.1570)), n = 4000, sd =  (0.1493))
plot(ts.sim, type = "l")
mag.sim <- arima.sim(mag.arima, 10000)


plot(mag.sub, type = "l")
fw <- fitdist(mag.hi, "weibull")

nobs <- length(mag.hi)
plot(sort(mag.hi), qweibull(1:nobs/nobs, fw$estimate[1], fw$estimate[2]))

ratvec <- u.hi / v.hi
plot(ratvec, type = "l")
thetavec <- atan(u.hi / v.hi)
plot(thetavec, type = "l")

theta.arima <- auto.arima(thetavec)
theta.resids <- theta.arima$residuals

theta.sim <- arima.sim(theta.arima, 10000)

plot(u.hi, type = "l")
plot(v.hi, type = "l")

# 
# plot(u.low, type = "l")
# plot(v.low, type = "l")
# plot(w.low, type = "l")
# plot(u.mid, type = "l")
# plot(v.mid, type = "l")
# plot(w.mid, type = "l")
# plot(u.hi, type = "l")
# plot(v.hi, type = "l")
# plot(w.hi, type = "l")


hist(u.hi)

n.obs <- length(u.hi)
dvec <- (1:n.obs)/(n.obs + 1)
qn <- qnorm(dvec)
#u.hi.arima <-auto.arima(u.hi)

df.seq <- 1:30

pval.seq <- rep(0, length(df.seq))
d.seq <- rep(0, length(df.seq))
mad.seq <- rep(0, length(df.seq))
log.like.seq <- rep(0, length(df.seq))
	
d.vec <- rep(0, length(df.seq))


u.hi.spline <- smooth.spline(u.hi, df = 30)
#u.hi.spline.resids <- u.hi - u.hi.spline$y
u.hi.spline.resids <- u.hi -mean(u.hi)
plot(u.hi, type = "l")
lines(u.hi.spline, col = "red")

mag.vec <- sqrt(u.hi**2 + v.hi**2 + w.hi**2)
dir.vec <- atan(u.hi / v.hi)

#plot(dir.vec[40000:50000], type = "l")

xtest.result <- read.table("xtestResult.txt", header = TRUE)
u.result <- xtest.result$u
v.result <- xtest.result$v

sd.u.hi <- sd(u.hi)
sd.v.hi <- sd(v.hi)
sd.result.u <- sd(u.result)
sd.result.v <- sd(v.result)

u.result.aug <- u.result * sd.u.hi / sd.result.u
v.result.aug <- v.result * sd.v.hi / sd.result.v

mean.u <- mean(u.hi)
mean.v <- mean(v.hi)

u.result <- u.result + mean.u
v.result <- v.result + mean.v

dir.result <- atan(u.result.aug / v.result.aug)
dir.uv <- atan(u.hi / v.hi)

uv.result <- data.frame(u = u.result.aug, v = v.result.aug)
cov.result <- cov(uv.result)
cov.result
obs.df <- data.frame(u = u.hi- mean(u.hi), v = v.hi - mean(v.hi))
obs.cov <- cov(obs.df)
obs.cov

inputcov <- matrix(c(3, -.5, -.5,  3), nrow  = 2)
input.eigen <- eigen(inputcov, symmetric = TRUE)

R <- t(input.eigen$vectors %*% (t(input.eigen$vectors) * sqrt(pmax(input.eigen$values, 0)))) ###dot prodoct of the eigenvectors *
###sqrt of the max of all eigenvalues or zero - whichever is greater

set.seed(12345)
X.init <- cbind(rt(1000, 30), rt(1000, 30))

X.mat <- as.matrix(X.init)

muvec <- colMeans(X.mat)
X2 <- X.mat - muvec
X <- X2
p <- length(muvec)

X1 <-  drop(muvec) +  input.eigen$vectors %*% diag(sqrt(pmax(input.eigen$values, 0)), p) %*% 
  t(X)
covX <- cov(Xt)


####U
for(cur.df in df.seq)
{
	sf <- 1/qt(.75, cur.df)
	temp.s <- mad(u.hi.spline.resids, sf)
	mad.seq[cur.df] <- temp.s
	med.temp <- median(temp.s)

	u.hi.nrm <- (u.hi.spline.resids - med.temp) / temp.s

	temp.lik <- dt(u.hi.nrm, df = cur.df)
	sum.lik <- sum(log(temp.lik))
	log.like.seq[cur.df] <- sum.lik
#   		
# 	fout <- paste("U-QQ-plot-", cur.df, ".png")
# 	png(fout)
# 		
# 	plot(qt(dvec, cur.df), sort(as.numeric(u.hi.nrm )))
# 
# 	dev.off()
	print(cur.df)
	
	u.scaled <- (u.hi.spline.resids - median(u.hi.spline.resids)) / mad(u.hi.spline.resids)
	
	ktest <- ks.test(u.scaled, "pt", cur.df)
	d.vec[cur.df] <- ktest$statistic
}

d.vec
which.min(d.vec)
par(mfrow = c(1, 1))

h1 <- hist(u.hi, xlim = c(-15, 15), breaks = 80)
b1 <- h1$breaks
####V


v.hi.spline <- smooth.spline(v.hi, df = 30)
#v.hi.spline.resids <- v.hi - v.hi.spline$y
v.hi.spline.resids <- v.hi - mean(v.hi)

mad.seq.v <- rep(NA, length(df.seq))


plot(v.hi, type = "l")
lines(v.hi.spline, col = "red")
abline(h = mean(v.hi), col = "red")
for(cur.df in df.seq)
{
  sf <- 1/qt(.75, cur.df)
  temp.s.v <- mad(v.hi.spline.resids, sf)
  mad.seq.v[cur.df] <- temp.s
  med.temp <- median(temp.s)
  
  v.hi.nrm <- (v.hi.spline.resids - med.temp) / temp.s
  
  temp.lik <- dt(v.hi.nrm, df = cur.df)
  sum.lik <- sum(log(temp.lik))
  log.like.seq.v[cur.df] <- sum.lik
  #   		
  # 	fout <- paste("U-QQ-plot-", cur.df, ".png")
  # 	png(fout)
  # 		
  # 	plot(qt(dvec, cur.df), sort(as.numeric(u.hi.nrm )))
  # 
  # 	dev.off()
  print(cur.df)
  
  v.scaled <- (v.hi.spline.resids - median(v.hi.spline.resids)) / mad(v.hi.spline.resids)
  
  ktest <- ks.test(v.scaled, "pt", cur.df)
  d.vec.v[cur.df] <- ktest$statistic
}


plot(d.vec.v)
opt.df.v <- which.min(d.vec.v)

u.hi.spline.diffinv <- diffinv(u.hi.spline.resids)
plot(u.hi.spline.diffinv, type = "l")
arima.u <- auto.arima(u.hi.spline.resids)

#arima.u <- auto.arima(u.hi.spline.resids, max.d = 0)

arima.v <- auto.arima(v.hi.spline.resids)

uvcov <- cov(u.hi.spline.resids, v.hi.spline.resids)
uvcor <- cor(u.hi.spline.resids, v.hi.spline.resids)
uvcor.test <- cor.test(u.hi.spline.resids, v.hi.spline.resids)
uvcor.test
uvcovmat <- cbind(c(1, uvcov), c(uvcov, 1))

residmat <- data.frame(ur = u.hi.spline.resids, vr = v.hi.spline.resids)
cov(residmat)


ar.residmat <- data.frame(ur = arima.u$residuals, v = arima.v$residuals)
cov(ar.residmat)
uvcor.test <- cor.test(arima.u$residuals, arima.v$residuals)

nobs.sim <- 1000

udiff <- diff(u.hi.spline.resids)
vdiff <- diff(v.hi.spline.resids)
usim <- arima.sim(arima.u, nobs.sim)

udiff.arima <- auto.arima(udiff)
vdiff.arima <- auto.arima(vdiff)


usim <- arima.sim(n = 1000, list(ar = c(0.5246  , -0.0731, 0.0788  ), ma = c( -0.3091   , -0.5881), sd = 2.398122))
usim <- arima.sim(n = 1000, udiff.arima)


vsim <- arima.sim(n = 1000, vdiff.arima)


###doing arima simulation directly

usim <- arima.sim(n = 1000, list(ar = c(.5  , .25), ma = c( -0.3   , -0.6)))
# 
# x1 <- rnorm(1000, 5, 3)
# x2 <- rnorm(1000, 2, 3)

x.init <- read.table("X-init.txt", header=  TRUE)
x1  <- x.init[,1]
x2 <- x.init[,2]
x12.cv <- matrix(c(3, -.5, -.5, 3), nrow = 2)


ev <- eigen(x12.cv, symmetric = TRUE) ###Get the eigenvalue/eigenvectors of the covariance matrix

R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0)))) ###dot prodoct of the eigenvectors *
###sqrt of the max of all eigenvalues or zero - whichever is greater
x.init <- read.table("X-init.txt", header=  TRUE) 
x1  <- x.init[,1]
x2 <- x.init[,2]
X <- matrix(cbind(x1, x2), ncol = 2)
head(X)
X.scale <- scale(X)
head(X.scale)
muvec <- colMeans(X)
sdvec <- apply(X, 2, "sd")
muvec
sdvec
X2 <- X - muvec
X <- X2
head(X)

p <- length(muvec)
X1 <-  drop(muvec) +  ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), p) %*% 
  t(X.scale)

X1 <- t(X1)
head(X1)
cov(X1)

x1.filt <- filter(X1[,1], c(1, c( -0.3   , -0.6)), sides = 1L)
x1.filt[seq_along( c(1, c( -0.3   , -0.6)))] <- 0
head(x1.filt)
x1.filt  <- filter(x1.filt, c(.5  , .25), method = "recursive")
head(x1.filt)

x2.filt <- filter(X1[,2], c(1, c( -0.3   , -0.6)), sides = 1L)
head(x2.filt)
x2.filt[seq_along( c(1, c( -0.3   , -0.6)))] <- 0
head(x2.filt)
x2.filt  <- filter(x2.filt, c(.5  , .25), method = "recursive")
head(x2.filt)

xmat.filt <- cbind(x1.filt, x2.filt)
#write.table(xmat.filt, "x-filt.txt", row.names = F)

