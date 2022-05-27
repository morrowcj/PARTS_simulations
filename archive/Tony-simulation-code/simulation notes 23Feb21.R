###############################
#implementing spatial autocofrrelation with panel regression method
###############################

n <- 4
p <- 5

w <- matrix(1,n,n) - diag(n)
W <- kronecker(diag(p), w)

lam <- .5
D <- solve(diag(n*p) - lam*W)
V <- D %*% t(D)
v <- diag(diag(V)^-.5)
(v %*% V %*% v)[1:6,1:6]

# queen
nSpace <- 21
m <- matrix(0, nSpace^2, nSpace^2)
index <- cbind(rep(1:nSpace, times = nSpace), rep(1:nSpace, each = nSpace))
for(i in 1:nSpace^2) {
	m[i,] <- (abs(index[i,1] - index[,1]) <=1) & (abs(index[i,2] - index[,2]) <= 1)
}
W <- m-diag(diag(m))
lam <- .12
D <- solve(diag(nSpace^2) - lam*W)
V <- D %*% t(D)
v <- diag(diag(V)^-.5)
#(v %*% V %*% v)

matrix((v %*% V %*% v)[ceiling(nSpace^2/2),], nrow=nSpace)
sum(matrix((v %*% V %*% v)[ceiling(nSpace^2/2),], nrow=nSpace))

plot(1:nSpace - ceiling(nSpace/2), matrix((v %*% V %*% v)[ceiling(nSpace^2/2),], nrow=nSpace)[ceiling(nSpace/2),])
curve(exp(-x/4), from=0, to=nSpace, add=T)


###############################
###############################
# simulation of a mixed model
###############################
###############################

library(mvtnorm)
library(Matrix)
library(geosphere)

#setwd("D:beautydata/arives/RSE ms/")

source("pglmm_ARMA_source_26Jan21.R")
source('remote_sensing_tools_28Dec20.R')

###############################
# simulation function
###############################

sim_ARMA21 <- function(mm, p, q, b, a, range, nugget, nSpace, cc, cc.sd = 0){
	xdim <- nSpace
	ydim <- nSpace
	nn <- nSpace^2

	location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
	colnames(location) <- c("lng", "lat")
	Dist <- distm(location)
	Dist <- Dist/(max(Dist)/2^.5)
	rownames(Dist) <- 1:nn
	colnames(Dist) <- 1:nn

	Vs <- nugget * diag(nrow(Dist)) + (1 - nugget) * exp(-Dist/range)

	m.start <- 100

	cc <- cc * matrix(1,nn,1) + rnorm(n=nn, sd = cc.sd)
	U <- rep(1:mm, times=nn)
	UU <- matrix(U/mean(U), ncol=mm, byrow=T)

	X <- matrix(NA, nrow=nn, ncol=mm+m.start)
	x <- matrix(0,nrow=nn, ncol=1)
	x1 <- matrix(0,nrow=nn, ncol=1)
	e1 <- t(rmvnorm(n=1, sigma=Vs))
	for(t in 1:(m.start + mm)){
		e <- t(rmvnorm(n=1, sigma=Vs))
		x.save <- x
		if(t <= m.start){
			x <- b1 * x +  b2 * x1 + e + a1*e1
		}else{
			x <- b1 * (x - cc*UU[,t - m.start]) +  b2 * (x1 - cc*UU[,t - m.start])+ cc*UU[,t - m.start] + e + a1*e1
		}
		x1 <- x.save
		e1 <- e
		X[,t] <- x
	}
	X <- X[,-(1:m.start)]

	dat <- data.frame(x=array(t(X)), time=rep(UU[1,], times=dim(X)[1]), site=rep(1:dim(X)[1], each=dim(X)[2]))
	dat$site <- as.factor(dat$site)
	return(list(dat = dat, Dist = Dist))
}

###############################
# examples
###############################

mm <- 30
p <- 2
q <- 1
b1 <- .5
b2 <- 0
a1 <- 0

nSpace <- 6
n <- nSpace^2
range <- .5
nugget <- .1

cc <- .2
cc.sd <- .5

sim <- sim_ARMA21(mm = mm, p = 2, q = 1, b = c(b1,b2), a = a1, range = range, nugget = nugget, nSpace = nSpace, cc = cc, cc.sd = cc.sd)
dat <- sim$dat
Dist <- sim$Dist
dat$site <- as.factor(dat$site)

#### Treating site as a random effect (but without spatial autocorrrelation
mod <- pglmm_ARMA(x ~ 1 + time, data = dat, p = 1, q = 0, distance.matrix = Dist, spatial.autocorrelation.function = "exponential", b = b1, a = 0, range = range, nugget = nugget, verbose = F, optimizer = "bobyqa", REML = F)
mod

#### Treating site as a random effect (but without spatial autocorrrelation
mod.re <- pglmm_ARMA(x ~ 1 + time + (1 | site) + (time | site), data = dat, p = 1, q = 0, distance.matrix = Dist, spatial.autocorrelation.function = "exponential", b = b1, range = range, nugget = nugget, verbose = F, optimizer = "bobyqa", REML = F)
mod.re


###############################
# PARTS
###############################
nugget.fit.flag <- T

X <- matrix(dat$x, byrow=T, ncol=mm)
t.scale <- dat$time[1:mm]

# fit AR-REML
dat.map <- AR_REML_fit(X,t.scale)
r.est <- spatialcor.fit(X, t.scale, Dist, fit.n.sample=nSpace^2, r.start=range, plot.fig=T, FUN = "exponential")

# construct the GLS correlation matrix and fit the GLS
V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")

if(nugget.fit.flag & r.est$spatialcor > .001){
	opt.nugget <- nugget.fit(formula='c ~ 1', dat.map, V, nugget.tol = 0.001)
	V <- (1-opt.nugget)*V + opt.nugget*diag(n)
}else{
	opt.nugget <- NA
	V <- V
}
if(r.est$spatialcor > .001){
	invcholV <- t(backsolve(chol(V), diag(n)))
} else {
	invcholV <- diag(n)
}

z <- GLS.fit(c ~ 1, data=dat.map, invcholV=invcholV)


mod
mod.re
c(z$coef, z$se, z$p.t)



###############################
###############################
# Simulation Study
###############################
###############################

FUN <- function(cc.sd, mm, b1, nSpace, range, nugget, cc, nrep = 1000){
	w <- data.frame(cc = cc, cc.sd = cc.sd)
	for(i.rep in 1:nrep){
		b2 <- 0
		a1 <- 0
		n <- nSpace^2


		sim <- sim_ARMA21(mm = mm, p = 2, q = 1, b = c(b1,b2), a = a1, range = range, nugget = nugget, nSpace = nSpace, cc = cc, cc.sd = cc.sd)
		dat <- sim$dat
		Dist <- sim$Dist
		dat$site <- as.factor(dat$site)

		#### Treating site as a random effect (but without spatial autocorrrelation
		mod <- pglmm_ARMA(x ~ 1 + time, data = dat, p = 1, q = 0, distance.matrix = Dist, spatial.autocorrelation.function = "exponential", b = b1, a = 0, range = range, nugget = nugget, verbose = F, optimizer = "bobyqa", REML = F)


		#### Treating site as a random effect (but without spatial autocorrrelation
		modre <- pglmm_ARMA(x ~ 1 + time + (1 | site) + (time | site), data = dat, p = 1, q = 0, distance.matrix = Dist, spatial.autocorrelation.function = "exponential", b = b1, range = range, nugget = nugget, verbose = F, optimizer = "bobyqa", REML = F)

		##### PARTS
		X <- matrix(dat$x, byrow=T, ncol=mm)
		t.scale <- dat$time[1:mm]

		# fit AR-REML
		dat.map <- AR_REML_fit(X,t.scale)
		r.est <- spatialcor.fit(X, t.scale, Dist, fit.n.sample=nSpace^2, r.start=range, plot.fig=F, FUN = "exponential")

		# construct the GLS correlation matrix and fit the GLS
		V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")

		# if(nugget.fit.flag & r.est$spatialcor > .001){
			# opt.nugget <- nugget.fit(formula='c ~ 1', dat.map, V, nugget.tol = 0.001)
			# V <- (1-opt.nugget)*V + opt.nugget*diag(n)
		# }else{
			# opt.nugget <- NA
			# V <- V
		# }
		# use input nugget
		V <- (1-nugget)*V + nugget*diag(n)

		if(r.est$spatialcor > .001){
			invcholV <- t(backsolve(chol(V), diag(n)))
		} else {
			invcholV <- diag(n)
		}

		z <- GLS.fit(c ~ 1, data=dat.map, invcholV=invcholV)

		# collect results
		w$mod.AIC <- mod$AIC
		w$mod.b <- mod$b
		w$mod.nugget <- mod$nugget
		w$mod.range <- mod$range
		w$mod.c <- mod$B[2]
		w$mod.c.zscore <- mod$B.zscore[2]
		w$mod.c.p <- mod$B.pvalue[2]

		w$modre.AIC <- modre$AIC
		w$modre.b <- modre$b
		w$modre.nugget <- modre$nugget
		w$modre.range <- modre$range
		w$modre.c <- modre$B[2]
		w$modre.c.zscore <- modre$B.zscore[2]
		w$modre.c.p <- modre$B.pvalue[2]
		w$modre.site <- modre$s2r[1]^.5
		w$modre.time.site <- modre$s2r[2]^.5

		w$gls.c <- z$coef
		w$gls.c.tscore <- z$t
		w$gls.c.p <- z$p.t
		w$gls.nugget <- nugget
		w$gls.range <- r.est$spatialcor

		if(i.rep == 1){
			write.table(w, file=paste0("pglmm_PARTS_b1=",b1,"_range=",range,"_nugget=",nugget,"_cc=",cc,"_cc.sd=",cc.sd,".csv"), sep=",", row.names=F, append=F)
		}else{
			write.table(w, file=paste0("pglmm_PARTS_b1=",b1,"_range=",range,"_nugget=",nugget,"_cc=",cc,"_cc.sd=",cc.sd,".csv"), sep=",", row.names=F, col.names = F, append=T)
		}
		#show(w)
		#show(modre)
	}
}
######################

mm <- 30
p <- 2
q <- 1
b1 <- .5
b2 <- 0
a1 <- 0

nSpace <- 6
n <- nSpace^2
range <- .3
nugget <- .1

cc <- .5
cc.sd <- .5

sim <- sim_ARMA21(mm = mm, p = 2, q = 1, b = c(b1,b2), a = a1, range = range, nugget = nugget, nSpace = nSpace, cc = cc, cc.sd = cc.sd)
dat <- sim$dat
X <- matrix(dat$x, byrow=T, ncol=mm)
matplot(t(X), lty=1, type="l")


#FUN(mm=mm, b1=b1, nSpace=nSpace, range=range, nugget=nugget, cc=cc, cc.sd=cc.sd, nrep = 10)


library(doParallel)

cc.sd.list <- .1*(0:10)
registerDoParallel(cl = 11)
getDoParWorkers()
foreach(cc.sd = cc.sd.list, .verbose=T, .packages=c("mgcv","geosphere","mvtnorm","Matrix")) %dopar% FUN(cc.sd, mm=mm, b1=b1, nSpace=nSpace, range=range, nugget=nugget, cc=cc, nrep = 500)
registerDoSEQ()



###############################
###############################
# Simulation Study for GLS with fitted nugget
###############################
###############################

FUN_GLS <- function(cc.sd, mm, b1, nSpace, range, nugget, nugget.gls, cc, nrep = 1000){
	w <- data.frame(cc = cc, cc.sd = cc.sd)
	for(i.rep in 1:nrep){
		b2 <- 0
		a1 <- 0
		n <- nSpace^2


		sim <- sim_ARMA21(mm = mm, p = 2, q = 1, b = c(b1,b2), a = a1, range = range, nugget = nugget, nSpace = nSpace, cc = cc, cc.sd = cc.sd)
		dat <- sim$dat
		Dist <- sim$Dist
		dat$site <- as.factor(dat$site)

		# fit AR-REML
		X <- matrix(dat$x, byrow=T, ncol=mm)
		t.scale <- dat$time[1:mm]

		dat.map <- AR_REML_fit(X,t.scale)
		r.est <- spatialcor.fit(X, t.scale, Dist, fit.n.sample=nSpace^2, r.start=range, plot.fig=F, FUN = "exponential")

		# construct the GLS correlation matrix and fit the GLS
		V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")

		if(is.null(nugget.gls) & r.est$spatialcor > .001){
			opt.nugget <- nugget.fit(formula='c ~ 1', dat.map, V, nugget.tol = 0.001)
		}else{
			opt.nugget <- nugget.gls
		}
		if(r.est$spatialcor > .001) V <- (1-opt.nugget)*V + opt.nugget*diag(n)

		if(r.est$spatialcor > .001){
			invcholV <- t(backsolve(chol(V), diag(n)))
		} else {
			invcholV <- diag(n)
		}

		z <- GLS.fit(c ~ 1, data=dat.map, invcholV=invcholV)

		# collect results

		w$gls.c <- z$coef
		w$gls.c.tscore <- z$t
		w$gls.c.p <- z$p.t
		w$gls.nugget <- opt.nugget
		w$gls.range <- r.est$spatialcor

		if(i.rep == 1){
			write.table(w, file=paste0("pglmm_PARTS_gls_b1=",b1,"_range=",range,"_nugget=",nugget,"_cc=",cc,"_cc.sd=",cc.sd,".csv"), sep=",", row.names=F, append=F)
		}else{
			write.table(w, file=paste0("pglmm_PARTS_gls_b1=",b1,"_range=",range,"_nugget=",nugget,"_cc=",cc,"_cc.sd=",cc.sd,".csv"), sep=",", row.names=F, col.names = F, append=T)
		}
		#show(w)
		#show(modre)
	}
}
######################

mm <- 30
p <- 2
q <- 1
b1 <- .5
b2 <- 0
a1 <- 0

nSpace <- 6
n <- nSpace^2
range <- .3
nugget <- .1
nugget.gls <- NULL

cc <- 0

library(doParallel)

cc.sd.list <- .1*(0:10)
registerDoParallel(cl = 6)
getDoParWorkers()
foreach(cc.sd = cc.sd.list, .verbose=T, .packages=c("mgcv","geosphere","mvtnorm","Matrix")) %dopar% FUN_GLS(cc.sd, mm=mm, b1=b1, nSpace=nSpace, range=range, nugget=nugget, nugget.gls=NULL, cc=cc, nrep = 500)
registerDoSEQ()

###############################
###############################
# Upload Simulation Study
###############################
###############################
p <- 2
q <- 1
b1 <- .5
b2 <- 0
a1 <- 0

nSpace <- 6
n <- nSpace^2
range <- .3
nugget <- .1

cc <- 0
cc.sd.list <- .1*(0:10)

alpha <- 0.05
alpha.gls <- .005
df <- data.frame(cc.sd = cc.sd.list)
count <- 0
for (cc.sd in cc.sd.list){
	count <- count + 1

	w <- read.csv(file=paste0("pglmm_ARMA_vs_PARTS_26Jan21/pglmm_PARTS_b1=",b1,"_range=",range,"_nugget=",nugget,"_cc=",cc,"_cc.sd=",cc.sd,".csv"))
	w1 <- read.csv(file=paste0("pglmm_ARMA_vs_PARTS_26Jan21/pglmm_PARTS_gls_b1=",b1,"_range=",range,"_nugget=",nugget,"_cc=",cc,"_cc.sd=",cc.sd,".csv"))

	df$mod.c[count] <- mean(w$mod.c.p < alpha)
	df$modre.c[count] <- mean(w$modre.c.p < alpha)
	df$gls.c[count] <- mean(w$gls.c.p < alpha)
	df$gls1.c[count] <- mean(w1$gls.c.p < alpha.gls)

	df$mod.dif.AIC[count] <- mean(w$modre.AIC - w$mod.AIC)

	df$time.site[count] <- mean(w$modre.time.site)

	df$mod.range[count] <- mean(w$mod.range)
	df$modre.range[count] <- mean(w$modre.range)
	df$gls.range[count] <- mean(w$gls.range)
	df$gls1.range[count] <- mean(w1$gls.range)

	df$mod.nugget[count] <- mean(w$mod.nugget)
	df$modre.nugget[count] <- mean(w$modre.nugget)
	df$gls.nugget[count] <- mean(w$gls.nugget)
	df$gls1.nugget[count] <- mean(w1$gls.nugget)

	df$mod.b[count] <- mean(w$mod.b)
	df$modre.b[count] <- mean(w$modre.b)
}
df0 <- df
df5 <- df

dft <- rbind(df0,df5)
dft
write.table(dft, file="pglmm_vs_PARTS.csv", sep=",", row.names=F, col.names = T)




