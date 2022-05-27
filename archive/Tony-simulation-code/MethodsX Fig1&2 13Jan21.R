# This file contains a lot of work to try to get a good method
# setwd("D:beautydata/arives/RSE ms/")

library(Matrix)
library(geosphere)
library(colorspace)

source('remote_sensing_tools_13Jan21.R')


###############################################################
# simulation for type I error and power
###############################################################

FUN <- function(bb, dat, do.GLS){
	for(r in r.list) {

		if(r > 0) {
			Dr <- as(t(chol(exp(-Dist/r))), "dgTMatrix")
		} else {
			Dr <- diag(nrow(Dist))
		}
	
		for(rep in 1:nrep){
			#set.seed(rep)
			
			output <- data.frame(r=r)
			output.p <- data.frame(r=r)
					
			e <- Dr %*% rnorm(n, sd = sd)
			b0 <- rnorm(n=n, 0, sd=b0.sd)
			y <- as.array(b0 + bb*dat$lat + e)
			dat$b0 <- b0
			dat$bb <- bb
			dat$y <- y
					
			if(do.GLS){
				# construct the GLS correlation matrix and fit the GLS
				V <- V.fit(Dist, spatialcor=r, FUN="exponential")
				
				if(nugget.fit.flag & r > .001){
					opt.nugget <- nugget.fit(formula='y ~ lat', dat, V, nugget.tol = 0.001)
					V <- (1-opt.nugget)*V + opt.nugget*diag(n)		
				}else{
					opt.nugget <- NA
					V <- V
				}
				if(r > .001){
					invcholV <- t(backsolve(chol(V), diag(n)))
				} else {
					invcholV <- diag(n)
				}
						
				z <- GLS.fit(y ~ lat, data=dat, invcholV=invcholV)
				
				output$F <- z$F
				output$p.F <- z$p.F
			}else{
				output$F <- NA
				output$p.F <- NA
			}
	
			# using partition
			z.p <- GLS.partition.data(y ~ 1 + lat, formula0 = y ~ 1, data=dat, spatialcor = r * r.correction, npart=npart, est.nugget=F, nugget.tol = 10^-3, nugget.interval = c(0,.5), spatial.autocor.FUN = "exponential")
			pvalue <- GLS.partition.pvalue(z.p, doFtest = T, nboot = 2000)
		
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
			output$rSSE <- z.p$rSSE
			output$rcoef <- z.p$rcoef[2]
	
			output$p.meanF <- pvalue$p.Fmean$pvalue
			output$p.meanF.rank <- pvalue$p.Fmean$rank.MSR
			output$p.chisq <- pvalue$p.chisq
			output$p.Fmin <- min(z.p$p.F.part)
			output$p.Fhochberg <- pvalue$p.Fhochberg
			output$p.fdr <- pvalue$p.Ffdr
			output$p.F1 <- z.p$p.F.part[1]
			output$p.F2 <- z.p$p.F.part[2]
			output$p.F3 <- z.p$p.F.part[3]
			
			output$p.t <- pvalue$p.t[2,4]
	
			# show(output)
	
			if(rep == 1){
				write.table(output, file=paste0("./X Fig1 results/X Fig1 r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F)
	
				write.table(output.rSSR.part, file=paste0("./X Fig1 results/X Fig1 rSSR.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.rSSE.part, file=paste0("./X Fig1 results/X Fig1 rSSE.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.SSR.part, file=paste0("./X Fig1 results/X Fig1 SSR.part r=", r, " n=", n, " bb", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.SSE.part, file=paste0("./X Fig1 results/X Fig1 SSE.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F)
				write.table(output.F.part, file=paste0("./X Fig1 results/X Fig1 F.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F)
			}else{
				write.table(output, file=paste0("./X Fig1 results/X Fig1 r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F, append=T)
	
				write.table(output.rSSR.part, file=paste0("./X Fig1 results/X Fig1 rSSR.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.rSSE.part, file=paste0("./X Fig1 results/X Fig1 rSSE.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.F.part, file=paste0("./X Fig1 results/X Fig1 F.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.SSR.part, file=paste0("./X Fig1 results/X Fig1 SSR.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F, append=T)
				write.table(output.SSE.part, file=paste0("./X Fig1 results/X Fig1 SSE.part r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv"), sep=",", row.names=F, col.names=F, append=T)
			}	
		}
	}
}

##############################
library(doParallel)


nSpace <- 40
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
max.untrans <- max(Dist)
Dist <- Dist/(max(Dist)/2^.5)
max.trans <- max(Dist)
c(max.untrans, max.trans)
r.correction <- max.untrans/max.trans/10^3

dat <- data.frame(lng=matrix(location[,1], ncol=1), lat=matrix(location[,2], ncol=1))

b0 <- 0
sd <- 1
b0.sd <- 0

fit.n.sample <- 400
nugget.fit.flag <- F
npart <- 8

r.list <- c(.03, .1)
# for nSpace = 60
#bb.list <- 2*(0:10)

# for nSpace = 40
bb.list <- 4*(0:10)

nrep <- 500

registerDoParallel(cl = 6)
getDoParWorkers()
foreach(bb = bb.list, .verbose=T, .packages=c("mgcv","geosphere")) %dopar% FUN(bb, dat=dat, do.GLS=TRUE)
registerDoSEQ()


##############################
# check figure
r <- .1
bb <- 0
Dr <- as(t(chol(exp(-Dist/r))), "dgTMatrix")
					
e <- Dr %*% rnorm(n, sd = sd)
b0 <- rnorm(n=n, 0, sd=b0.sd)
y <- as.array(b0 + bb*dat$lat + e)
dat$b0 <- b0
dat$bb <- bb
dat$y <- y

image(matrix(dat$y, nSpace, nSpace))
z <- GLS.partition.data(y ~ 1 + lat, formula0 = y ~ 1, data=dat, spatialcor = r*r.correction, npart=npart, est.nugget=F, nugget.tol = 10^-3, nugget.interval = c(0,.5), spatial.autocor.FUN = "exponential")
GLS.partition.pvalue(z, doFtest = T, nboot = 2000)
	

##############################
# summarize
r.list <- c(.03,.1)

b0.sd <- 0

n <- 3600
npart <- 8
bb.list <- 2*(0:10)
# n <- 1600
# npart <- 8
# bb.list <- 4*(0:10)

alpha <- .05

d <- data.frame(r=rep(NA, length(r.list)*length(cc.list)), n=n)
count <- 0
for(r in r.list) for(bb in bb.list) {
	count <- count + 1
	d$r[count] <- r
	d$bb[count] <- bb
	filename <- paste0("X Fig1 r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv")
	if(is.element(filename, dir("./X Fig1 results"))){
		filename <- paste0("./X Fig1 results/X Fig1 r=", r, " n=", n, " bb=", bb, " b0.sd=", b0.sd, " npart=", npart, " 10Jan21.csv")
		output <- read.csv(filename)
		d$r.est[count] <- mean(output$r)	
		d$reject.F[count] <- mean(output$p.F < alpha)
		d$reject.p.F1[count] <- mean(output$p.F1 < alpha)
		#d$reject.p.F2[count] <- mean(output$p.F2 < alpha)
		#d$reject.p.F3[count] <- mean(output$p.F3 < alpha)
		d$reject.p.Fmin[count] <- mean(output$p.Fmin < alpha)
		d$reject.p.meanF[count] <- mean(output$p.meanF < alpha)
		d$reject.p.chisq[count] <- mean(output$p.chisq < alpha)
		d$reject.p.Fhochberg[count] <- mean(output$p.Fhochberg < alpha)
		d$reject.p.fdr[count] <- mean(output$p.fdr < alpha)
		d$reject.p.t[count] <- mean(output$p.t < alpha)
	}
}

d.05 <- d
d.01 <- d


# alpa = 0.05
# d40
      # r    n bb r.est reject.F reject.p.F1 reject.p.Fmin reject.p.meanF reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.t
# 1  0.03 1600  0  0.03    0.042       0.052         0.248          0.042          0.044              0.032        0.036      0.038
# 2  0.03 1600  4  0.03    0.116       0.092         0.366          0.126          0.122              0.082        0.096      0.130
# 3  0.03 1600  8  0.03    0.294       0.184         0.540          0.296          0.302              0.230        0.242      0.306
# 4  0.03 1600 12  0.03    0.560       0.314         0.758          0.518          0.528              0.368        0.398      0.522
# 5  0.03 1600 16  0.03    0.818       0.540         0.916          0.796          0.798              0.674        0.710      0.802
# 6  0.03 1600 20  0.03    0.922       0.674         0.974          0.918          0.918              0.814        0.850      0.922
# 7  0.03 1600 24  0.03    0.994       0.842         0.994          0.992          0.990              0.944        0.954      0.992
# 8  0.03 1600 28  0.03    0.996       0.918         0.998          0.992          0.994              0.976        0.982      0.994
# 9  0.03 1600 32  0.03    1.000       0.982         1.000          1.000          1.000              0.996        0.998      1.000
# 10 0.03 1600 36  0.03    1.000       0.998         1.000          1.000          1.000              1.000        1.000      1.000
# 11 0.03 1600 40  0.03    1.000       0.998         1.000          1.000          1.000              1.000        1.000      1.000
# 12 0.10 1600  0  0.10    0.064       0.070         0.170          0.064          0.068              0.026        0.040      0.070
# 13 0.10 1600  4  0.10    0.060       0.056         0.146          0.058          0.056              0.022        0.032      0.056
# 14 0.10 1600  8  0.10    0.098       0.102         0.178          0.098          0.098              0.048        0.058      0.096
# 15 0.10 1600 12  0.10    0.124       0.112         0.220          0.126          0.128              0.050        0.064      0.130
# 16 0.10 1600 16  0.10    0.210       0.196         0.364          0.204          0.214              0.128        0.148      0.204
# 17 0.10 1600 20  0.10    0.308       0.266         0.432          0.270          0.280              0.180        0.208      0.278
# 18 0.10 1600 24  0.10    0.434       0.368         0.562          0.398          0.404              0.258        0.330      0.396
# 19 0.10 1600 28  0.10    0.568       0.470         0.708          0.526          0.538              0.396        0.454      0.528
# 20 0.10 1600 32  0.10    0.628       0.566         0.762          0.610          0.614              0.462        0.508      0.614
# 21 0.10 1600 36  0.10    0.784       0.654         0.842          0.734          0.732              0.620        0.648      0.728
# 22 0.10 1600 40  0.10    0.842       0.770         0.914          0.802          0.800              0.692        0.746      0.798

# d60
      # r    n bb r.est   reject.F reject.p.F1 reject.p.Fmin reject.p.meanF reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.t
# 1  0.03 3600  0  0.03 0.04582485  0.05193483     0.1904277     0.04684318     0.04684318         0.03462322   0.03971487 0.04786151
# 2  0.03 3600  2  0.03 0.08146640  0.08452138     0.2556008     0.07535642     0.07637475         0.04378819   0.05193483 0.07942974
# 3  0.03 3600  4  0.03 0.19653768  0.14663951     0.3951120     0.18737271     0.18839104         0.13340122   0.15071283 0.18533605
# 4  0.03 3600  6  0.03 0.36456212  0.28207739     0.5600815     0.35234216     0.34928717         0.25763747   0.28716904 0.35030550
# 5  0.03 3600  8  0.03 0.58553971  0.44297352     0.7627291     0.56008147     0.57026477         0.43380855   0.47148676 0.56211813
# 6  0.03 3600 10  0.03 0.79633401  0.64154786     0.9083503     0.77800407     0.78004073         0.67107943   0.70468432 0.78411405
# 7  0.03 3600 12  0.03 0.90529532  0.78309572     0.9541752     0.88594705     0.88900204         0.80855397   0.83401222 0.89002037
# 8  0.03 3600 14  0.03 0.96435845  0.86558045     0.9887984     0.95621181     0.95824847         0.91242363   0.92973523 0.95824847
# 9  0.03 3600 16  0.03 0.98778004  0.94501018     0.9949084     0.98065173     0.98370672         0.96639511   0.97352342 0.98167006
# 10 0.03 3600 18  0.03 0.99592668  0.97963340     0.9969450     0.99490835     0.99490835         0.98574338   0.98879837 0.99490835
# 11 0.03 3600 20  0.03 1.00000000  0.99592668     1.0000000     1.00000000     1.00000000         0.99796334   1.00000000 1.00000000
# 12 0.10 3600  0  0.10 0.05511022  0.06112224     0.1172345     0.05511022     0.05711423         0.02004008   0.03306613 0.05911824
# 13 0.10 3600  2  0.10 0.05505506  0.05505506     0.1151151     0.05305305     0.05205205         0.02502503   0.03203203 0.05305305
# 14 0.10 3600  4  0.10 0.07815631  0.07414830     0.1503006     0.07715431     0.08116232         0.03807615   0.05010020 0.08116232
# 15 0.10 3600  6  0.10 0.10610611  0.10110110     0.1841842     0.10010010     0.09809810         0.05705706   0.07407407 0.10010010
# 16 0.10 3600  8  0.10 0.14829659  0.13026052     0.2274549     0.14228457     0.14028056         0.07414830   0.09018036 0.14028056
# 17 0.10 3600 10  0.10 0.19739479  0.18136273     0.3006012     0.18937876     0.18837675         0.10821643   0.13426854 0.18837675
# 18 0.10 3600 12  0.10 0.29058116  0.25450902     0.3917836     0.27655311     0.27655311         0.17535070   0.21142285 0.27655311
# 19 0.10 3600 14  0.10 0.34454638  0.30886850     0.4668705     0.31192661     0.32212029         0.20591233   0.25076453 0.31906218
# 20 0.10 3600 16  0.10 0.41182467  0.36799185     0.5249745     0.38634047     0.38430173         0.27115189   0.30988787 0.38430173
# 21 0.10 3600 18  0.10 0.51580020  0.47298675     0.6279307     0.48521916     0.48725790         0.33129460   0.38226300 0.48521916
# 22 0.10 3600 20  0.10 0.60652396  0.54841998     0.7145770     0.59021407     0.60142712         0.43628950   0.50050968 0.59734964

# alpa = 0.01
# d40
      # r    n bb r.est reject.F reject.p.F1 reject.p.Fmin reject.p.meanF reject.p.chisq reject.p.Fhochberg reject.p.fdr reject.p.t
# 1  0.03 1600  0  0.03    0.010       0.014         0.060          0.006          0.006              0.004        0.004      0.006
# 2  0.03 1600  4  0.03    0.034       0.008         0.116          0.038          0.038              0.018        0.020      0.038
# 3  0.03 1600  8  0.03    0.130       0.076         0.274          0.126          0.128              0.100        0.102      0.126
# 4  0.03 1600 12  0.03    0.314       0.132         0.456          0.282          0.290              0.162        0.182      0.292
# 5  0.03 1600 16  0.03    0.626       0.278         0.736          0.580          0.580              0.426        0.450      0.580
# 6  0.03 1600 20  0.03    0.822       0.466         0.872          0.782          0.798              0.608        0.636      0.790
# 7  0.03 1600 24  0.03    0.954       0.656         0.966          0.928          0.930              0.812        0.834      0.934
# 8  0.03 1600 28  0.03    0.978       0.782         0.986          0.974          0.976              0.936        0.944      0.976
# 9  0.03 1600 32  0.03    1.000       0.918         1.000          0.998          0.998              0.982        0.990      0.998
# 10 0.03 1600 36  0.03    1.000       0.968         1.000          1.000          1.000              1.000        1.000      1.000
# 11 0.03 1600 40  0.03    1.000       0.988         1.000          1.000          1.000              1.000        1.000      1.000
# 12 0.10 1600  0  0.10    0.018       0.006         0.040          0.010          0.012              0.004        0.008      0.014
# 13 0.10 1600  4  0.10    0.018       0.010         0.036          0.016          0.016              0.012        0.014      0.014
# 14 0.10 1600  8  0.10    0.018       0.030         0.060          0.024          0.026              0.006        0.010      0.026
# 15 0.10 1600 12  0.10    0.030       0.030         0.070          0.030          0.030              0.018        0.022      0.026
# 16 0.10 1600 16  0.10    0.076       0.060         0.150          0.056          0.074              0.044        0.050      0.072
# 17 0.10 1600 20  0.10    0.130       0.110         0.210          0.122          0.128              0.074        0.084      0.126
# 18 0.10 1600 24  0.10    0.196       0.154         0.314          0.190          0.196              0.126        0.148      0.196
# 19 0.10 1600 28  0.10    0.314       0.252         0.446          0.290          0.302              0.192        0.226      0.294
# 20 0.10 1600 32  0.10    0.406       0.330         0.526          0.370          0.374              0.300        0.320      0.374
# 21 0.10 1600 36  0.10    0.542       0.440         0.658          0.492          0.516              0.406        0.436      0.498
# 22 0.10 1600 40  0.10    0.642       0.522         0.740          0.608          0.618              0.470        0.516      0.620

# d60
      # r    n bb r.est    reject.F reject.p.F1 reject.p.Fmin reject.p.meanF reject.p.chisq reject.p.Fhochberg reject.p.fdr  reject.p.t
# 1  0.03 3600  0  0.03 0.010183299 0.012219959    0.04887984    0.009164969    0.011201629        0.010183299  0.010183299 0.012219959
# 2  0.03 3600  2  0.03 0.021384929 0.014256619    0.06211813    0.020366599    0.019348269        0.014256619  0.017311609 0.019348269
# 3  0.03 3600  4  0.03 0.072301426 0.040733198    0.16293279    0.065173116    0.064154786        0.048879837  0.052953157 0.062118126
# 4  0.03 3600  6  0.03 0.173116090 0.110997963    0.30957230    0.167006110    0.169042770        0.121181263  0.135437882 0.167006110
# 5  0.03 3600  8  0.03 0.343177189 0.217922607    0.50814664    0.327902240    0.325865580        0.243380855  0.258655804 0.318737271
# 6  0.03 3600 10  0.03 0.585539715 0.373727088    0.73523422    0.550916497    0.557026477        0.435845214  0.473523422 0.553971487
# 7  0.03 3600 12  0.03 0.746435845 0.550916497    0.83706721    0.707739308    0.711812627        0.605906314  0.628309572 0.708757637
# 8  0.03 3600 14  0.03 0.889002037 0.703665988    0.93584521    0.856415479    0.862525458        0.771894094  0.791242363 0.858452138
# 9  0.03 3600 16  0.03 0.962321792 0.841140530    0.97657841    0.948065173    0.951120163        0.897148676  0.917515275 0.952138493
# 10 0.03 3600 18  0.03 0.981670061 0.913441955    0.99083503    0.977596741    0.979633401        0.960285132  0.965376782 0.979633401
# 11 0.03 3600 20  0.03 0.997963340 0.968431772    1.00000000    0.997963340    0.996945010        0.991853360  0.992871690 0.997963340
# 12 0.10 3600  0  0.10 0.009018036 0.009018036    0.02304609    0.006012024    0.008016032        0.002004008  0.002004008 0.008016032
# 13 0.10 3600  2  0.10 0.014014014 0.013013013    0.02402402    0.013013013    0.013013013        0.006006006  0.010010010 0.014014014
# 14 0.10 3600  4  0.10 0.022044088 0.021042084    0.04609218    0.022044088    0.022044088        0.012024048  0.013026052 0.023046092
# 15 0.10 3600  6  0.10 0.031031031 0.026026026    0.06306306    0.030030030    0.034034034        0.013013013  0.017017017 0.033033033
# 16 0.10 3600  8  0.10 0.044088176 0.044088176    0.08817635    0.042084168    0.043086172        0.018036072  0.028056112 0.042084168
# 17 0.10 3600 10  0.10 0.070140281 0.069138277    0.12124248    0.064128257    0.066132265        0.042084168  0.051102204 0.065130261
# 18 0.10 3600 12  0.10 0.122244489 0.091182365    0.19839679    0.119238477    0.117234469        0.072144289  0.082164329 0.116232465
# 19 0.10 3600 14  0.10 0.152905199 0.140672783    0.23139653    0.146788991    0.148827727        0.096839959  0.113149847 0.150866463
# 20 0.10 3600 16  0.10 0.199796126 0.168195719    0.28644241    0.194699286    0.193679918        0.117227319  0.145769623 0.192660550
# 21 0.10 3600 18  0.10 0.254841998 0.218144750    0.36085627    0.243628950    0.241590214        0.160040775  0.188583078 0.236493374
# 22 0.10 3600 20  0.10 0.367991845 0.304791030    0.47298675    0.331294597    0.339449541        0.235474006  0.266055046 0.337410805

######################################################
# X Fig 1

pdf(paste0("X Fig1 npart=",npart," n=",n," alpha = .05 and .01 10Jan21.pdf"), height=8, width=8)

	par(mfrow=c(2,2), mai=c(.8,.8,.5,.2))
	
	xlab=rep(expression(beta), 2)
	ylab=c("Proportion Rejected","")
	mlab=c("a", "b", "c", "d")
	mainlab=c(expression(italic(r) ~ '= 0.03,' ~ alpha ~ '= 0.05'),expression(italic(r) ~ '= 0.1,' ~ alpha ~ '= 0.05'),expression(italic(r) ~ '= 0.03,' ~ alpha ~ '= 0.01'),expression(italic(r) ~ '= 0.1,' ~ alpha ~ '= 0.01'))
	col.list <- c("black","red","red","red", "orange", "orange","blue")
	lty.list <- c(1,1,2,3,1,2,1)
	
	##########
	alpha = 0.05
	d <- d.05
	r <- r.list[1]
	i <- 0
	for(j in c(5,8,9,12,10,11,6)){
		i <- i+1
		dd <- d[d$r==r,]
		if(i==1) {
			plot(dd[,j] ~ dd$bb, typ="l", ylim=c(0,1), xlab=xlab[2], ylab=ylab[1], cex.lab=1.2, col=col.list[i])
			mtext(mlab[1], side=4, at=1, adj=-.5, las=1, cex=1.2)
			mtext(mainlab[1], side=3, at=10, las=1, cex=1.2)
		} else {
			lines(dd[,j] ~ bb, data=dd, col=col.list[i], lty=lty.list[i])
		}
	}
	lines(c(0,100), c(alpha,alpha), lty=3)
	
	r <- r.list[2]
	i <- 0
	for(j in c(5,8,9,12,10,11,6)){
		i <- i+1
		dd <- d[d$r==r,]
		if(i==1) {
			plot(dd[,j] ~ dd$bb, typ="l", ylim=c(0,1), xlab=xlab[1], ylab=ylab[1], cex.lab=1.2, col=col.list[i])
			mtext(mlab[2], side=4, at=1, adj=-.5, las=1, cex=1.2)
			mtext(mainlab[2], side=3, at=10, las=1, cex=1.2)
		} else {
			lines(dd[,j] ~ bb, data=dd, col=col.list[i], lty=lty.list[i])
		}
	}
	lines(c(0,100), c(alpha,alpha), lty=3)
	
	##########
	alpha = 0.01
	d <- d.01
	r <- r.list[1]
	i <- 0
	for(j in c(5,8,9,12,10,11,6)){
		i <- i+1
		dd <- d[d$r==r,]
		if(i==1) {
			plot(dd[,j] ~ dd$bb, typ="l", ylim=c(0,1), xlab=xlab[2], ylab=ylab[1], cex.lab=1.2, col=col.list[i])
			mtext(mlab[3], side=4, at=1, adj=-.5, las=1, cex=1.2)
			mtext(mainlab[3], side=3, at=10, las=1, cex=1.2)
		} else {
			lines(dd[,j] ~ bb, data=dd, col=col.list[i], lty=lty.list[i])
		}
	}
	lines(c(0,100), c(alpha,alpha), lty=3)
	
	r <- r.list[2]
	i <- 0
	for(j in c(5,8,9,12,10,11,6)){
		i <- i+1
		dd <- d[d$r==r,]
		if(i==1) {
			plot(dd[,j] ~ dd$bb, typ="l", ylim=c(0,1), xlab=xlab[1], ylab=ylab[1], cex.lab=1.2, col=col.list[i])
			mtext(mlab[4], side=4, at=1, adj=-.5, las=1, cex=1.2)
			mtext(mainlab[4], side=3, at=10, las=1, cex=1.2)
		} else {
			lines(dd[,j] ~ bb, data=dd, col=col.list[i], lty=lty.list[i])
		}
	}
	lines(c(0,100), c(alpha,alpha), lty=3)


	if(is.na(dd$reject.p.F1[1])){
		legend(x=.0, y=1, legend=c(expression(italic(P[F])),expression(italic(P[LRT])),expression(italic(P[t])),expression(italic(P[hoch])), expression(italic(P[fdr])), expression(italic(P[single]))), col=col.list[-1], lty=lty.list[-1], box.col="white")
	}else{
		legend(x=.0, y=1, legend=c(expression(italic(P[GLS])),expression(italic(P[F])),expression(italic(P[LRT])),expression(italic(P[t])),expression(italic(P[hoch])), expression(italic(P[fdr])), expression(italic(P[single]))), col=col.list, lty=lty.list, box.col="white")
	}

dev.off()


######################################################
######################################################
# Fig. 2
######################################################
######################################################

nSpace <- 60
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
max.untrans <- max(Dist)
Dist <- Dist/(max(Dist)/2^.5)
max.trans <- max(Dist)
c(max.untrans, max.trans)
r.correction <- max.untrans/max.trans/10^3

dat <- data.frame(lng=matrix(location[,1], ncol=1), lat=matrix(location[,2], ncol=1))

b0 <- 0
sd <- 1
b0.sd <- 0

fit.n.sample <- 400
nugget.fit.flag <- F
npart <- 8

r <- .03
#r <- .1
# for r=0.03
bb <- 5
# for r=0.1
#bb <- 20
Dr <- as(t(chol(exp(-Dist/r))), "dgTMatrix")
					
e <- Dr %*% rnorm(n, sd = sd)
b0 <- rnorm(n=n, 0, sd=b0.sd)
y <- as.array(b0 + bb*dat$lat + e)
dat$b0 <- b0
dat$bb <- bb
dat$y <- y

nrep <- 100
output <- data.frame(rep=1:nrep)
for(i in 1:nrep){
	z <- GLS.partition.data(y ~ 1 + lat, formula0 = y ~ 1, data=dat, spatialcor = r*r.correction, npart=npart, est.nugget=F, nugget.tol = 10^-3, nugget.interval = c(0,.5), spatial.autocor.FUN = "exponential", min.num.cross.part = npart)
	
	output$p[i] <- correlated.chisq(Fmean.obs = z$Fmean, rSSR = z$rSSR, df1 = z$df1, npart = z$npart)$pvalue
	
	show(output[i,])
}

hist(output$p)
output003 <- output
output01 <- output

write.table(output003, file="MethodsX Fig2 0.03.csv", sep=",", row.names=F)
write.table(output01, file="MethodsX Fig2 0.1.csv", sep=",", row.names=F)


pdf(paste0("X Fig2 npart=",npart," n=",n," 13Jan21.pdf"), height=4, width=8)

	par(mfrow=c(1,2), mai=c(.8,.8,.5,.2))
	
	mlab=c("a", "b", "c", "d")
	mainlab=c(expression(italic(r) ~ '= 0.03'),expression(italic(r) ~ '= 0.1'))
	
	hist(output003$p, breaks=.002*(0:40)+.02, main='', ylab="Frequency", xlab=expression(italic(P[LRT])), freq=T)
	
	mtext(mlab[1], side=4, at=13, adj=-.5, las=1, cex=1.3)
	mtext(mainlab[1], side=3, at=.08, las=1, cex=1.2, padj=5)
	
	hist(output01$p, breaks=.002*(0:40)+.02, main='', ylab="Frequency", xlab=expression(italic(P[LRT])), freq=T)
	
	mtext(mlab[2], side=4, at=18, adj=-.5, las=1, cex=1.3)
	mtext(mainlab[2], side=3, at=.08, las=1, cex=1.2, padj=5)
dev.off()


######################################################
######################################################
# investigate full vs. average rSSR
######################################################
######################################################

nSpace <- 40
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
max.untrans <- max(Dist)
Dist <- Dist/(max(Dist)/2^.5)
max.trans <- max(Dist)
c(max.untrans, max.trans)
r.correction <- max.untrans/max.trans/10^3

dat <- data.frame(lng=matrix(location[,1], ncol=1), lat=matrix(location[,2], ncol=1))

b0 <- 0
sd <- 1
b0.sd <- 0

fit.n.sample <- 400
nugget.fit.flag <- F
npart <- 8

r <- .1
bb <- 40
Dr <- as(t(chol(exp(-Dist/r))), "dgTMatrix")
					
e <- Dr %*% rnorm(n, sd = sd)
b0 <- rnorm(n=n, 0, sd=b0.sd)
y <- as.array(b0 + bb*dat$lat + e)
dat$b0 <- b0
dat$bb <- bb
dat$y <- y

source('remote_sensing_tools_13Jan21.R')

nrep <- 10
output <- data.frame(rep=1:nrep)
for(i in 1:nrep){
	z <- GLS.partition.data(y ~ 1 + lat, formula0 = y ~ 1, data=dat, spatialcor = r*r.correction, npart=npart, est.nugget=F, nugget.tol = 10^-3, nugget.interval = c(0,.5), spatial.autocor.FUN = "exponential", min.num.cross.part = npart)
	
	output$rSSR[i] <- correlated.chisq(Fmean.obs = z$Fmean, rSSR = z$rSSR, df1 = z$df1, npart = z$npart)$pvalue
	output$rSSR.part[i] <- correlated.chisq(Fmean.obs = z$Fmean, rSSR = z$rSSR.part, df1 = z$df1, npart = z$npart)$pvalue
}
output

   # rep       rSSR  rSSR.part
# 1    1 0.01794429 0.01794249
# 2    2 0.01571913 0.01571599
# 3    3 0.01312957 0.01312823
# 4    4 0.01497800 0.01497637
# 5    5 0.01338951 0.01338548
# 6    6 0.01634019 0.01633848
# 7    7 0.02204575 0.02204432
# 8    8 0.01579339 0.01579179
# 9    9 0.02051542 0.02051366
# 10  10 0.01805002 0.01804728




