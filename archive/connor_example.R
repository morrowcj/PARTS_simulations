library(ggplot2)
set.seed(8)

# 1) Data setup ----

# paramters for simulation
## time steps
ntime = 9
## map size
map.cols = map.rows = 20
## total number of pixels
npix = map.cols * map.rows
## spatial correlation range parameter
r = .25 # 25% the width of the map
## nugget
nug = .10
## diagonal nugget matrix
nugmat = diag(rep(nug, npix))
## pixel coordinates
coords = scales::rescale(as.matrix(expand.grid(x = as.double(1:map.cols),
                                               y = as.double(1:map.rows))), to = c(0, 1))
## Distance between pixels (scaled between 0 and 1)
D <- scales::rescale(as.matrix(dist(coords)), to = c(0, 1))

# 2) Covariance ----
## We use the exponential covariance function:
covar_func <- function(d, range = r){
  covar = exp(-d/range)
  return(covar)
}

## Semivariogram:
  # For the true covariance matrix, we also add in a nugget that absorbs
  # additional variation. A traditional way to visualize spatial variance is
  # with a semivariogram. The semivariance which is simply 1 - covariance. Here's a plot of the
  # function:
curve(1 - ((1-nug) * covar_func(x)), from = 0, to = 1,
      ylab = "semivariance", xlab = "distance", main = "Semivariogram",
      ylim = c(0,1))
# add nugget line and legend
abline(h = nug, col = "forestgreen", lty = 2)
legend(x = .75, y = .5, lty = 2, col = "forestgreen", legend = "nugget")

## Here's the full covariance matrix for our data
Sigma = nugmat + (1 - nug) * covar_func(D, r)

## a Spatial error structure is created by multiplying Sigma by a random variable
del = scale(Sigma %*% rnorm(npix))

## visualize the spatial error on a map
ggplot(data.frame(coords, del), aes(x = x, y = y, fill = del)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "x coordinate", y = "y coordinate", fill = "value")


## 3) Estimate r from the covariance matrix: ----

## NOTE The next step is estimating the range parameter.
  # Theoretically, the correlation among CLS residuals model should be close to Sigma.
  # So, in remote parts, this estimation is done with that residual correlation matrix (an estimate of Sigma).
  # but, for simplicity, we use the true covariance matrix here.

## First, we need to combine the covars and distances into a data frame
v.cov.resids = Sigma[upper.tri(Sigma, diag = TRUE)] # vector of true covariance
W = data.frame(dist = v.dist, cor = v.cov.resids) # paired dist-cov values

## Then estimate r with ML
r.est = coef(fit <- nls(cor ~ exp(-dist/r), data = W, start = list(r = .1), trace = FALSE))
r.est # very close to true value

## We also estimate the nugget similarly, though I won't show it here.
nug.est = .139 # arbitrary nugget estimate

# So, our estimate for spatial covariance that we would use in GLS (excluding the nugget) would be
Sigma.est = diag(rep(nug.est, npix)) + (1 - nug.est) * covar_func(D, r.est)

## We can see the average difference between our estimate and the true Sigma is small
mean(Sigma - Sigma.est)
