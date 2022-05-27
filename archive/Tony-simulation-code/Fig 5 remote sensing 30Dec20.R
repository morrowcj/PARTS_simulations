#################################################
#################################################
# see MethodsX Fig1 16Dec20.R
#################################################
#################################################

#setwd("D:beautydata/arives/RSE ms/")

library(doParallel)
library(Matrix)
library(geosphere)
library(colorspace)

source('remote_sensing_tools_28Dec20.R')

###############################################################
# simulation for type I error and power
###############################################################

FUN <- function(cc, dat, location, nrep, r.list, npart.list, do.GLS){
	for(r in r.list) for(npart in npart.list){

		if(r > 0) {
			Dr <- as(t(chol(exp(-Dist/r))), "dgTMatrix")
		} else {
			Dr <- diag(nrow(Dist))
		}
	
		for(rep in 1:nrep){
			
			output <- data.frame(r=r)
					
			dat$c0 <- rnorm(n=n, 0, sd=c0.sd)
			X <- simX('~ 0 + c0 + landclass', data=dat, coef=c(1, cc * slopes), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=rep)
				
			# fit AR-REML
			dat.map <- AR_REML_fit(X,t.scale)
			
			# this adds the landscape variable to dat.map
			dat.map$landclass <- dat$landclass			
			dat.map$lng <- location[,1]
			dat.map$lat <- location[,2]
	
			r.est <- spatialcor.fit(X, t.scale, Dist, fit.n.sample=fit.n.sample, r.start=r, plot.fig=F, FUN = "exponential")
			output$r <- r.est$spatialcor
			output$r.sigma <- r.est$spatialcor.sigma
					
			time1 <- Sys.time()
			
			if(do.GLS){
				# construct the GLS correlation matrix and fit the GLS
				V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")
				
				if(nugget.fit.flag & r.est$spatialcor > .001){
					opt.nugget <- nugget.fit(formula='c ~ landclass', dat.map, V, nugget.tol = 0.001)
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
						
				z <- GLS.fit(c ~ 0 + landclass, data=dat.map, invcholV=invcholV)
				z0 <- GLS.fit(c ~ 1, data=dat.map, invcholV=invcholV)
				
				names.coef <- names(z$coef)
				output.coef <- matrix(z$coef, nrow=1)
				colnames(output.coef) <- names.coef
				
				output$GLSp.F <- z$p.F				
				output$GLSp.t0 <- z0$p.t
				output$GLSnugget <- opt.nugget
			}else{
				output$GLSp.F <- NA			
				output$GLSp.t0 <- NA
				output$GLSnugget <- NA
			}
	
			time2 <- Sys.time()
			output$time1 <- as.numeric(time2-time1)
					
			# using partition			
			r.est <- spatialcor.fit.data(X, t.scale, data=dat.map, fit.n.sample=fit.n.sample, r.start=r, plot.fig=F, FUN = "exponential")

			z.p <- GLS.partition.data(c ~ landclass, formula0 = c ~ 1, data=dat.map, spatialcor = r.est$spatialcor, npart=npart, est.nugget=T, nugget.tol = 10^-3, nugget.interval = c(0,.5), spatial.autocor.FUN = "exponential", min.num.cross.part = 5)
			pvalue <- GLS.partition.pvalue(z.p, doFtest = F, nboot=10^4)
			z.p0 <- GLS.partition.data(c ~ 1, formula0 = c ~ 0, data=dat.map, spatialcor = r.est$spatialcor, npart=npart, est.nugget=T, nugget.tol = 10^-3, nugget.interval = c(0,.5), spatial.autocor.FUN = "exponential", min.num.cross.part = 5)
			pvalue0 <- GLS.partition.pvalue(z.p0, doFtest = F)
		
			names.coef <- names(z.p$coef)
			output.coef.part <- matrix(z.p$coef, nrow=1)
			colnames(output.coef.part) <- names.coef
			
			output.SSR.part <- matrix(z.p$SSR.part, nrow=1)
			output.SSE.part <- matrix(z.p$SSE.part, nrow=1)
			output.F.part <- matrix(z.p$F.part, nrow=1)
			output.rSSR.part <- matrix(z.p$rSSR.part[!is.na(z.p$rSSR.part)], nrow=1)
			output.rSSE.part <- matrix(z.p$rSSE.part[!is.na(z.p$rSSE.part)], nrow=1)
			output$nugget.p <- z.p$nugget
			
			output$rSSR <- z.p$rSSR
			output$rSSR <- z.p$rSSR
			output$rSSE <- z.p$rSSE
			output$rcoef <- z.p$rcoef[2]
	
			output$df1 <- z.p$df1
			output$p.chisq <- pvalue$p.chisq
			output$rank.deficient.MSR <- pvalue$rank.deficient.MSR
			output$p.Fmean <- pvalue$p.Fmean[1]$pvalue
			
			output$p.Fmin <- min(z.p$p.F.part)
			output$p.Fhochberg <- pvalue$p.Fhochberg
			output$p.fdr <- pvalue$p.Ffdr
			output$p.F1 <- z.p$p.F.part[1]
			output$p.F2 <- z.p$p.F.part[2]
			output$p.F3 <- z.p$p.F.part[3]
			output$p.t <- pvalue$p.t[2,4]
	
			output$r.part <- r.est$spatialcor
			
			output$df10 <- z.p0$df1
			output$p.chisq0 <- pvalue0$p.chisq
			output$p.Fmin0 <- min(z.p0$p.F.part)
			output$p.Fhochberg0 <- pvalue0$p.Fhochberg
			output$p.fdr0 <- pvalue0$p.Ffdr
			output$p.F10 <- z.p0$p.F.part[1]
			output$p.F20 <- z.p0$p.F.part[2]
			output$p.F30 <- z.p0$p.F.part[3]
			output$p.t0 <- pvalue0$p.t[4]

			time3 <- Sys.time()
			output$time2 <- as.numeric(time3-time2)
			
			#show(output)
			if(rep == 1){
				write.table(output, file=paste0("Fig5 r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F)
	
				write.table(output.coef.part, file=paste0("Fig5 output.coef.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.rSSR.part, file=paste0("Fig5 rSSR.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.rSSE.part, file=paste0("Fig5 rSSE.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.SSR.part, file=paste0("Fig5 SSR.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.SSE.part, file=paste0("Fig5 SSE.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.F.part, file=paste0("Fig5 F.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F)
			}else{
				write.table(output, file=paste0("Fig5 r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
	
				write.table(output.coef.part, file=paste0("Fig5 output.coef.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.rSSR.part, file=paste0("Fig5 rSSR.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.rSSE.part, file=paste0("Fig5 rSSE.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.F.part, file=paste0("Fig5 F.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.SSR.part, file=paste0("Fig5 SSR.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.SSE.part, file=paste0("Fig5 SSE.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv"), sep=",", row.names=F, col.names=F, append=T)
			}	
		}
	}
}

##############################

nSpace <- 40
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
image(landscape)
dat <- data.frame(landclass=matrix(landscape, ncol=1) - 1)
dat$landclass <- as.factor(dat$landclass)
dat$lat <- location[,2]

n.obs <- 30
n.burn <- 10

b0 <- 0
b <- .2
s <- 1
c0.sd <- .4

t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

slopes <- 0:3

fit.n.sample <- 400
nugget.fit.flag <- T

npart.list <- 8
r.list <- c(.03, .1)

# for nSpace = 40
cc.list <- .02*(0:10)

do.GLS <- TRUE

registerDoParallel(cl = 11)
getDoParWorkers()
foreach(cc = cc.list, .verbose=T, .packages=c("mgcv","geosphere")) %dopar% FUN(cc, dat=dat, location=location, nrep=nrep, r.list=r.list, npart.list=npart.list, do.GLS=do.GLS)
registerDoSEQ()



#################
# test type I errors
cc <- 0
do.GLS <- F
nrep <- 500

npart.list <- c(4,8,12)
registerDoParallel(cl = 3)
getDoParWorkers()
foreach(npart.list = npart.list, .verbose=T, .packages=c("mgcv","geosphere")) %dopar% FUN(cc, dat=dat, location=location, nrep=nrep, r.list=r.list, npart.list=npart.list, do.GLS=do.GLS)
registerDoSEQ()
#################




######################################################
# summarize
cc.list <- .02*(0:10)
r.list <- c(.03,.1)

c0.sd <- .4

n <- 1600
npart <- 8

alpha <- .05

d <- data.frame(r=rep(NA, length(r.list)*length(cc.list)), n=n)
count <- 0
for(r in r.list) for(cc in cc.list) {
	count <- count + 1
	d$r[count] <- r
	d$cc[count] <- cc
	filename <- paste0("Fig5 r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 30Dec20.csv")
	if(is.element(filename, dir("./Fig5 results"))){
		filename <- paste0("./Fig5 results/", filename)
		output <- read.csv(filename)
		d$r.est[count] <- mean(output$r)	
		d$reject.p.GLS[count] <- mean(output$GLSp.F < alpha)
		d$reject.p.chisq[count] <- mean(output$p.chisq < alpha)
		d$reject.p.Fhochberg[count] <- mean(output$p.Fhochberg < alpha)
		d$reject.p.fdr[count] <- mean(output$p.fdr < alpha)
		d$reject.p.F1[count] <- mean(output$p.F1 < alpha)
				
		d$reject.p.GLS0[count] <- mean(output$GLSp.t0 < alpha)
		d$reject.p.chisq0[count] <- mean(output$p.chisq0 < alpha)
		d$reject.p.Fhochberg0[count] <- mean(output$p.Fhochberg0 < alpha)
		d$reject.p.fdr0[count] <- mean(output$p.fdr0 < alpha)
		d$reject.p.t0[count] <- mean(output$p.t0 < alpha)
		d$reject.p.t0[count] <- mean(output$p.t0 < alpha)
		d$reject.p.F10[count] <- mean(output$p.F10 < alpha)
	}
}
d

      # r    n   cc      r.est reject.p.GLS reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.F1 reject.p.GLS0
# 1  0.03 1600 0.00 0.02909084        0.060          0.094              0.078        0.080       0.062         0.050
# 2  0.03 1600 0.02 0.02909084        0.076          0.118              0.090        0.096       0.062         0.090
# 3  0.03 1600 0.04 0.02909084        0.136          0.170              0.124        0.128       0.076         0.222
# 4  0.03 1600 0.06 0.02909084        0.248          0.268              0.166        0.170       0.114         0.394
# 5  0.03 1600 0.08 0.02909084        0.406          0.394              0.242        0.262       0.178         0.630
# 6  0.03 1600 0.10 0.02909084        0.588          0.542              0.354        0.370       0.238         0.792
# 7  0.03 1600 0.12 0.02909084        0.778          0.714              0.504        0.530       0.332         0.910
# 8  0.03 1600 0.14 0.02909084        0.904          0.826              0.650        0.678       0.424         0.982
# 9  0.03 1600 0.16 0.02909084        0.962          0.918              0.780        0.798       0.514         0.990
# 10 0.03 1600 0.18 0.02909084        0.984          0.964              0.864        0.884       0.626         0.996
# 11 0.03 1600 0.20 0.02909084        0.998          0.984              0.932        0.946       0.730         1.000
# 12 0.10 1600 0.00 0.09677207        0.064          0.090              0.060        0.066       0.066         0.056
# 13 0.10 1600 0.02 0.09677207        0.076          0.088              0.066        0.068       0.070         0.058
# 14 0.10 1600 0.04 0.09677207        0.116          0.118              0.072        0.074       0.076         0.064
# 15 0.10 1600 0.06 0.09677207        0.204          0.158              0.088        0.094       0.106         0.098
# 16 0.10 1600 0.08 0.09677207        0.366          0.236              0.124        0.134       0.144         0.144
# 17 0.10 1600 0.10 0.09677207        0.556          0.324              0.188        0.208       0.208         0.186
# 18 0.10 1600 0.12 0.09677207        0.738          0.448              0.298        0.308       0.258         0.228
# 19 0.10 1600 0.14 0.09677207        0.850          0.568              0.390        0.418       0.336         0.294
# 20 0.10 1600 0.16 0.09677207        0.930          0.698              0.532        0.562       0.420         0.376
# 21 0.10 1600 0.18 0.09677207        0.972          0.776              0.662        0.684       0.522         0.462
# 22 0.10 1600 0.20 0.09677207        0.992          0.860              0.768        0.782       0.598         0.542
   # reject.p.chisq0 reject.p.Fhochberg0 reject.p.fdr0 reject.p.t0 reject.p.F10
# 1            0.062               0.050         0.052       0.062        0.052
# 2            0.092               0.078         0.084       0.082        0.076
# 3            0.214               0.148         0.158       0.216        0.128
# 4            0.418               0.260         0.276       0.418        0.202
# 5            0.616               0.444         0.474       0.626        0.346
# 6            0.774               0.632         0.674       0.796        0.484
# 7            0.902               0.808         0.826       0.898        0.636
# 8            0.966               0.908         0.920       0.974        0.750
# 9            0.988               0.964         0.974       0.988        0.850
# 10           0.994               0.982         0.984       0.996        0.934
# 11           1.000               0.990         0.992       1.000        0.968
# 12           0.066               0.028         0.038       0.064        0.064
# 13           0.064               0.032         0.038       0.062        0.064
# 14           0.078               0.044         0.054       0.074        0.070
# 15           0.096               0.064         0.076       0.092        0.086
# 16           0.132               0.084         0.096       0.126        0.112
# 17           0.176               0.108         0.128       0.172        0.160
# 18           0.224               0.142         0.176       0.214        0.202
# 19           0.280               0.186         0.212       0.272        0.254
# 20           0.354               0.238         0.268       0.342        0.314
# 21           0.432               0.294         0.348       0.430        0.392
# 22           0.504               0.366         0.426       0.496        0.456
######################################################
# check with larger map
######################################################

###############
# nSpace = 40
# npart=8
     # r    n cc      r.est reject.p.GLS reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.F1 reject.p.GLS0 reject.p.chisq0 reject.p.Fhochberg0 reject.p.fdr0 reject.p.t0
# 1 0.03 1600  0 0.02909084   0.05600000     0.08800000         0.06200000    0.0660000  0.05400000        0.0580          0.0680          0.04800000    0.05200000      0.0580
# 2 0.10 1600  0 0.09704491   0.07421875     0.08203125         0.05078125    0.0546875  0.08984375        0.0625          0.0625          0.03515625    0.04296875      0.0625

# npart = 4
     # r    n cc      r.est reject.p.GLS reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.F1 reject.p.GLS0 reject.p.chisq0 reject.p.Fhochberg0 reject.p.fdr0 reject.p.t0
# 1 0.03 1600  0 0.02909084        0.058          0.086              0.066        0.066       0.050         0.056           0.054               0.044         0.050       0.056
# 2 0.10 1600  0 0.09677207        0.048          0.090              0.070        0.074       0.052         0.054           0.058               0.036         0.042       0.054

###############
# nSpace = 60

# npart=12
     # r    n cc      r.est reject.p.GLS reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.F1 reject.p.GLS0 reject.p.chisq0 reject.p.Fhochberg0 reject.p.fdr0 reject.p.t0
# 1 0.03 3600  0 0.02905861           NA          0.054              0.068        0.074        0.06            NA           0.050               0.036          0.04       0.054
# 2 0.10 3600  0 0.09732488           NA          0.048              0.058        0.062        0.06            NA           0.058               0.022          0.03       0.054

# npart=8
     # r    n cc      r.est reject.p.GLS reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.F1 reject.p.GLS0 reject.p.chisq0 reject.p.Fhochberg0 reject.p.fdr0 reject.p.t0
# 1 0.03 3600  0 0.02905861           NA          0.052              0.046        0.046       0.044            NA           0.052               0.042         0.046       0.048
# 2 0.10 3600  0 0.09732488           NA          0.046              0.034        0.038       0.060            NA           0.054               0.022         0.030       0.052

# npart=4
     # r    n cc      r.est reject.p.GLS reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.F1 reject.p.GLS0 reject.p.chisq0 reject.p.Fhochberg0 reject.p.fdr0 reject.p.t0
# 1 0.03 3600  0 0.02905861           NA          0.044              0.036        0.038       0.036            NA           0.046               0.028         0.028       0.050
# 2 0.10 3600  0 0.09732488           NA          0.044              0.038        0.042       0.032            NA           0.046               0.032         0.036       0.046


######################################################
# Fig 5
######################################################

pdf(paste0("Fig5 remote sensing npart=",npart," n=",n," alpha =", alpha," 30Dec20.pdf"), height=7, width=8)

	par(mfcol=c(2,2), mai=c(.8,.8,.1,.4))
	
	xlab=rep(expression(zeta), 2)
	ylab=c("Proportion Rejected","")
	mlab=c("(a)", "(b)", "(c)", "(d)")
	mainlab=c(expression(italic(r) ~ '= 0.03'),expression(italic(r) ~ '= 0.1'))
	col.list <- c("black","red", "orange", "orange","blue")
	lty.list <- c(1,1,1,2,1)
	
	i.lab <- 0
	for (r in r.list) {
		i.lab <- i.lab+1
		i <- 0
		for(j in c(5,6,7,8,9)){
			i <- i+1
			dd <- d[d$r==r,]
			if(i==1) {
				plot(dd[,j] ~ dd$cc, typ="l", ylim=c(0,1), xlab=xlab[2], ylab=ylab[1], cex.lab=1.2, col=col.list[i])
				mtext(mlab[i.lab], side=4, at=1, adj=-.5, las=1, cex=1.2)
				mtext(mainlab[1], side=3, at=10, las=1, cex=1.2)
			} else {
				lines(dd[,j] ~ cc, data=dd, col=col.list[i], lty=lty.list[i])
			}
		}
		lines(c(0,100), c(alpha,alpha), lty=3)
	}
	
	i.lab <- 2
	for (r in r.list) {
		i.lab <- i.lab+1
		i <- 0
		for(j in c(10,11,12,13,15)){
			i <- i+1
			dd <- d[d$r==r,]
			if(i==1) {
				plot(dd[,j] ~ dd$cc, typ="l", ylim=c(0,1), xlab=xlab[1], ylab=ylab[1], cex.lab=1.2, col=col.list[i])
				mtext(mlab[i.lab], side=4, at=1, adj=-.5, las=1, cex=1.2)
				mtext(mainlab[2], side=3, at=10, las=1, cex=1.2)
			} else {
				lines(dd[,j] ~ cc, data=dd, col=col.list[i], lty=lty.list[i])
			}
		}
		lines(c(0,100), c(alpha,alpha), lty=3)
	}
	
	legend(x=.0, y=1, legend=c(expression(italic(P[GLS])),expression(italic(P[part])),expression(italic(P[hoch])), expression(italic(P[fdr])), expression(italic(P[single]))), col=col.list, lty=lty.list, box.col="white")

dev.off()


######################################################
# investigate corSSE and corSSR
r <- .03
cc <- 0
n <- 1600
npart <- 8
c0.sd <- .25
	
filename <- paste0("Fig5  r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 18Nov20.csv")
output <- read.csv(filename, header=T)

filename <- paste0("Fig5  SSR.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 18Nov20.csv")
SSR <- read.csv(filename, header=F)
filename <- paste0("Fig5  rSSR.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 18Nov20.csv")
rSSR <- read.csv(filename, header=F)
cor(SSR)
(sum(cor(SSR)) - npart)/(npart*(npart-1))
mean(as.matrix(rSSR))

filename <- paste0("Fig5  SSE.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 18Nov20.csv")
SSE <- read.csv(filename, header=F)
filename <- paste0("Fig5  rSSE.part r=", r, " n=", n, " cc=", cc, " c0.sd=", c0.sd, " npart=", npart, " 18Nov20.csv")
rSSE <- read.csv(filename, header=F)
cor(SSE)
(sum(cor(SSE)) - npart)/(npart*(npart-1))
mean(as.matrix(rSSE))


hist(as.numeric(as.matrix(SSR)))

