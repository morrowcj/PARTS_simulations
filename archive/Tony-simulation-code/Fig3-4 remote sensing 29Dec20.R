library(mvtnorm)
library(lattice)
library(raster)
library(Kendall)
library(Matrix)
library(MASS)
library(geosphere)
library(mblm)
library(fields)

library(colorspace)

source('remote_sensing_tools_28Dec20.R')

###############################
# Set up spatial covariance matrix
###############################

# set up variables
nSpace <- 100
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
Dist <- Dist/(max(Dist)/2^.5)

# distribution of 4 land classes
n.size <- nSpace/4
n.cluster <- 2
landscape <- kronecker(kronecker(matrix(1,ncol=n.cluster,nrow=n.cluster), matrix(1:4, nrow=2, ncol=2)), matrix(1, nrow=n.size, ncol=n.size))
dat <- data.frame(landclass=matrix(landscape, ncol=1) - 1)
dat$landclass <- as.factor(dat$landclass)

n.obs <- 30
n.burn <- 10

b0 <- 0
b <- .2
s <- 1
c0.sd <- 0

slopes <- 0:3

t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

nugget.fit.flag <- T

alpha <- .1

slopes <- 0:3
cc <- 0
dat$c0 <- rnorm(n=n, 0, sd=c0.sd)

#################
# Fig 3 AR-reml only
r.list <- c(0,.1,.1)

seed <- 1
for(i in 1:3){
	r <- r.list[i]
	if(i == 3) seed <- 2

	if(r > 0) {
		Dr <- t(chol(exp(-Dist/r)))
	} else {
		Dr <- diag(nrow(Dist))
	}
	X <- simX('~0 + c0', data=dat, coef=c(1), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=seed)
	
	w <- data.frame(site=1:n)
	for(counter in 1:dim(X)[1]){
		x <- X[counter,]
		w$r[counter] <- r
		w$b[counter] <- b
		w$c[counter] <- cc
		
		z <- AR_reml(x ~ t.scale)
		w$cc[counter] <- z$beta[2]
		w$p[counter] <- z$Pr[2]
		w$b[counter] <- z$b
	}
	if(i == 1) w0 <- w
	if(i == 2) w1 <- w
	if(i == 3) w2 <- w
}

pdf(file=paste0("Fig3 remote sensing b=", b, " cc=", cc, " n.obs=", n.obs, " alpha=", alpha, " 29Dec20.pdf"), height=4, width=12)
	par(mfrow = c(1,3), mai=c(.1,.1,.4,.55))
	col.pal <- hcl.colors(21, palette="Green-Brown", rev = TRUE)
	
	z <- w0

	show(mean(z$p < alpha))

	zmax <- max(abs(z$cc))
	
	matrix.fig <- matrix(z$cc, ncol=ydim)
	mask <- matrix(z$p, ncol=ydim)
	matrix.fig[mask > alpha] <- NA
	image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
	mtext(side = 3, expression(paste(italic(r), " = 0")), cex=1.5)
	mtext(side=3, "(a)", at = 1.05, cex=1.5)

	z <- w1

	show(mean(z$p < alpha))

	zmax <- max(abs(z$cc))
	
	matrix.fig <- matrix(z$cc, ncol=ydim)
	mask <- matrix(z$p, ncol=ydim)
	matrix.fig[mask > alpha] <- NA
	image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
	mtext(side=3, expression(paste(italic(r), " = 0.1")), cex=1.5)
	mtext(side=3, "(b)", at = 1.05, cex=1.5)
	
	z <- w2

	show(mean(z$p < alpha))

	zmax <- max(abs(z$cc))
	
	matrix.fig <- matrix(z$cc, ncol=ydim)
	mask <- matrix(z$p, ncol=ydim)
	matrix.fig[mask > alpha] <- NA
	image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
	mtext(side=3, expression(paste(italic(r), " = 0.1")), cex=1.5)
	mtext(side=3, "(c)", at = 1.05, cex=1.5)
dev.off()


###############################
#Fig4: land classes 
###############################


#################################################
# set up variables
nSpace <- 100
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
Dist <- Dist/(max(Dist)/2^.5)

# distribution of 4 land classes
n.size <- nSpace/4
n.cluster <- 2
landscape <- kronecker(kronecker(matrix(1,ncol=n.cluster,nrow=n.cluster), matrix(1:4, nrow=2, ncol=2)), matrix(1, nrow=n.size, ncol=n.size))
dat <- data.frame(landclass=matrix(landscape, ncol=1) - 1)
dat$landclass <- as.factor(dat$landclass)

n.obs <- 30
n.burn <- 10

b0 <- 0
b <- .2
s <- 1
c0.sd <- .2

slopes <- 0:3

t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

nugget.fit.flag <- T

r <- .1

if(r > 0) {
	Dr <- t(chol(exp(-Dist/r)))
} else {
	Dr <- diag(nrow(Dist))
}

##################################
# simulation

# for cc = 0	
seed <- 10

slopes <- 0:3
cc <- 0
dat$c0 <- rnorm(n=n, 0, sd=c0.sd)

X <- simX('~0 + c0 + landclass', data=dat, coef=c(1, cc * slopes), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=seed)

dat.map <- AR_REML_fit(X,t.scale)
dat.map$landscape <- as.factor(dat$landclass)
dat.map$lng <- location[,1]
dat.map$lat <- location[,2]

# look for correlation between c and b
par(mfrow=c(1,3), mai=c(.8,.8,.3,.3))
hist(dat.map$b)
lines(c(b,b),c(0,2000),col="red")
hist(dat.map$c)
plot(c ~ b, data=dat.map)

sd(dat.map$c)

fit.n.sample <- 2000
r.est <- spatialcor.fit(X, t.scale, Dist=Dist, fit.n.sample=fit.n.sample, plot.fig=T)
r.est

V <- exp(-Dist/r.est$spatialcor)
nugget <- nugget.fit(formula='c ~ 0 + landscape', dat.map, V, nugget.tol = 0.00001, verbose = T)
nugget

Vn <- (1 - nugget) * V + nugget * diag(n)

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z <- GLS.fit(c ~ 0 + landscape, data=dat.map, invcholV=invcholVn)

# save output
c.data <- aggregate(dat.map$c, by=list(landscape=dat.map$landscape), FUN=mean)[,2]
se.data <- aggregate(dat.map$c, by=list(landscape=dat.map$landscape), FUN=mean)[,2]/(n/4)^.5
matrix.fig <- matrix(dat.map$c, ncol=ydim)
mask <- matrix(dat.map$p.c, ncol=ydim)

M <- list(c.data=c.data, se.data=se.data, matrix.fig=matrix.fig, mask=mask, coef=z$coef, se=z$se, FF=z$F, p.FF=z$p.F, varcov=z$varcov, r=r.est$spatialcor, nugget=nugget, coef0=z$coef0, se0=z$se0)

saveRDS(file=paste0("Fig4 cc=",cc," seed=",seed," 29Dec19.RDS"), M)

# for cc = 0.1
seed <- 10

slopes <- 0:3
cc <- .1
dat$c0 <- rnorm(n=n, 0, sd=c0.sd)

X <- simX('~0 + c0 + landclass', data=dat, coef=c(1, cc * slopes), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=seed)

dat.map <- AR_REML_fit(X,t.scale)
dat.map$landscape <- as.factor(dat$landclass)
dat.map$lng <- location[,1]
dat.map$lat <- location[,2]

# look for correlation between c and b
par(mfrow=c(1,3), mai=c(.8,.8,.3,.3))
hist(dat.map$b)
lines(c(b,b),c(0,2000),col="red")
hist(dat.map$c)
plot(c ~ b, data=dat.map)

sd(dat.map$c)

fit.n.sample <- 2000
r.est <- spatialcor.fit(X, t.scale, Dist=Dist, fit.n.sample=fit.n.sample, plot.fig=T)
r.est

V <- exp(-Dist/r.est$spatialcor)
nugget <- nugget.fit(formula='c ~ 0 + landscape', dat.map, V, nugget.tol = 0.00001, verbose = T)
nugget

Vn <- (1 - nugget) * V + nugget * diag(n)

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z <- GLS.fit(c ~ 0 + landscape, data=dat.map, invcholV=invcholVn)

# save output
c.data <- aggregate(dat.map$c, by=list(landscape=dat.map$landscape), FUN=mean)[,2]
se.data <- aggregate(dat.map$c, by=list(landscape=dat.map$landscape), FUN=mean)[,2]/(n/4)^.5
matrix.fig <- matrix(dat.map$c, ncol=ydim)
mask <- matrix(dat.map$p.c, ncol=ydim)

M <- list(c.data=c.data, se.data=se.data, matrix.fig=matrix.fig, mask=mask, coef=z$coef, se=z$se, FF=z$F, p.FF=z$p.F, varcov=z$varcov, r=r.est$spatialcor, nugget=nugget, coef0=z$coef0, se0=z$se0)

saveRDS(file=paste0("Fig4 cc=",cc," seed=",seed," 29Dec19.RDS"), M)


##################################
# plotting
Fig4.plot <- function(M, lab){
	col.pal <- hcl.colors(21, palette="Green-Brown", rev = TRUE)
	land.pal <- hcl.colors(21, palette="Sunset",)[13:21]

	matrix.fig <- M$matrix.fig
	mask <- M$mask
	c.data <- M$c.data
	se.data <- M$se.data
	c.est <- M$c.est
	se.est <- M$se.est
	
	par(mai=c(.22,.22,.22,.35))
	
	matrix.fig[mask > alpha] <- NA
	landscape.fig <- landscape
	landscape.fig[mask < alpha] <- NA
	image(landscape.fig, xaxt="n", yaxt="n", col=land.pal)
	mtext(lab[1], side=4, at=.98, adj=-.5, las=1, cex=1.2)
	image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, add=T, zlim=c(-zmax,zmax))
	text("1", x=.125, y=.125, cex=1.5)
	text("2", x=.375, y=.125, cex=1.5)
	text("3", x=.125, y=.375, cex=1.5)
	text("4", x=.375, y=.375, cex=1.5)
	
	par(mai=c(.8,.8,.3,.35))
	
	arg <- 1:4
	plot(arg+.15, c.data, xlim=c(.5,4.5), ylim=c(-.4,.5), ylab="Time trend", xlab="Land class", main="", pch=1, cex.lab=1.5, col="red")
	mtext(lab[2], side=4, at=.5, adj=-.5, las=1, cex=1.2)
	arrows(x0=arg+.15, y0=c.data-max(se.data, .0001), y1=c.data+max(se.data, .0001), angle=90, code=3, length=.05, col="red")
	
	# coef
	points(arg+.05, M$coef, col="black", pch=15)
	arrows(x0=arg, y0=M$coef-M$se, y1=M$coef+M$se, angle=90, length=.05, code=3)	
	
	# for conditional confidence intervals
	S <- M$varcov
	p <- length(M$coef)
	se.cond <- array(p)
	for(i in 1:p){
		S11 <- S[i,i]
		S22 <- S[-i,-i]
		S12 <- S[i,-i]
		se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
	}
	arrows(x0=arg+.1, y0=M$coef - se.cond, y1=M$coef + se.cond, angle=90, code=3, length=.05, col="blue")

}

# set up variables
nSpace <- 100
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

n.size <- nSpace/4
n.cluster <- 2
landscape <- kronecker(kronecker(matrix(1,ncol=n.cluster,nrow=n.cluster), matrix(1:4, nrow=2, ncol=2)), matrix(1, nrow=n.size, ncol=n.size))

n.obs <- 30
n.burn <- 10

b0 <- 0
b <- .2
s <- 1
c0.sd <- .2

r <- .1
seed <- 10

cc <- 0
M1 <- readRDS(file="Fig4 cc=0 seed=10 29Dec19.RDS")

cc <- .1
M2 <- readRDS(file="Fig4 cc=0.1 seed=10 29Dec19.RDS")

pdf(file=paste0("Fig4 landscape r=",r, " cc=", cc, " b=", b, " n.obs=", n.obs, " alpha=", alpha, " seed=",seed," 29Dec19.pdf"), height=8, width=8)
par(mfcol=c(2,2))

zmax <- max(abs(M1$matrix.fig), abs(M2$matrix.fig))
alpha <- .1

Fig4.plot(M1, lab=c("(a)","(b)"))
Fig4.plot(M2, lab=c("(c)","(d)"))

dev.off()

M1$FF
M1$p.FF
M1$coef0/M1$se0
2*min(pt(abs(M1$coef0/M1$se0), df=10000, lower.tail=T),pt(abs(M1$coef0/M1$se0), df=10000, lower.tail=F))
M1$r
M1$nugget

anova(lm(matrix(M1$matrix.fig,ncol=1) ~ as.factor(matrix(landscape,ncol=1))))
summary(lm(matrix(M1$matrix.fig,ncol=1) ~ 1))


M2$FF
M2$p.FF
M2$coef0/M2$se0
2*min(pt(abs(M2$coef0/M2$se0), df=10000, lower.tail=T),pt(abs(M2$coef0/M2$se0), df=10000, lower.tail=F))
M2$r
M2$nugget

anova(lm(matrix(M2$matrix.fig,ncol=1) ~ as.factor(matrix(landscape,ncol=1))))
summary(lm(matrix(M2$matrix.fig,ncol=1) ~ 1))

# > M1$FF
# [1] 0.9563446
# > M1$p.FF
# [1] 0.4123046
# > M1$coef0/M1$se0
           # [,1]
# [1,] -0.1260782
# > 2*min(pt(abs(M1$coef0/M1$se0), df=10000, lower.tail=T),pt(abs(M1$coef0/M1$se0), df=10000, lower.tail=F))
# [1] 0.8996725
# > M1$r
        # r 
# 0.1034756 
# > M1$nugget
# [1] 0.05950413
# > 
# > anova(lm(matrix(M1$matrix.fig,ncol=1) ~ as.factor(matrix(landscape,ncol=1))))
# Analysis of Variance Table

# Response: matrix(M1$matrix.fig, ncol = 1)
                                         # Df Sum Sq Mean Sq F value    Pr(>F)    
# as.factor(matrix(landscape, ncol = 1))    3  332.3 110.780  159.86 < 2.2e-16 ***
# Residuals                              9996 6927.2   0.693                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > summary(lm(matrix(M1$matrix.fig,ncol=1) ~ 1))

# Call:
# lm(formula = matrix(M1$matrix.fig, ncol = 1) ~ 1)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -2.70593 -0.60611 -0.01455  0.63433  2.55301 

# Coefficients:
             # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.056025   0.008521  -6.575  5.1e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8521 on 9999 degrees of freedom

# > 
# > 
# > M2$FF
# [1] 12.68464
# > M2$p.FF
# [1] 2.849218e-08
# > M2$coef0/M2$se0
         # [,1]
# [1,] 0.743298
# > 2*min(pt(abs(M2$coef0/M2$se0), df=10000, lower.tail=T),pt(abs(M2$coef0/M2$se0), df=10000, lower.tail=F))
# [1] 0.4573187
# > M2$r
        # r 
# 0.1034756 
# > M2$nugget
# [1] 0.06063522
# > 
# > anova(lm(matrix(M2$matrix.fig,ncol=1) ~ as.factor(matrix(landscape,ncol=1))))
# Analysis of Variance Table

# Response: matrix(M2$matrix.fig, ncol = 1)
                                         # Df Sum Sq Mean Sq F value    Pr(>F)    
# as.factor(matrix(landscape, ncol = 1))    3  217.0  72.322   104.5 < 2.2e-16 ***
# Residuals                              9996 6917.8   0.692                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > summary(lm(matrix(M2$matrix.fig,ncol=1) ~ 1))

# Call:
# lm(formula = matrix(M2$matrix.fig, ncol = 1) ~ 1)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -2.87630 -0.60293 -0.00461  0.62951  2.64916 

# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.093084   0.008447   11.02   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8447 on 9999 degrees of freedom

