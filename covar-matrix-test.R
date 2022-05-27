map.width = 16 # number of columns (in pixels)
ntime = 5 # number of time points
nug = 0
r = .5
# pixels in the map
npix = map.width ^ 2
# coordinates of each pixel
coords = scales::rescale(as.matrix(expand.grid(x = as.double(1:map.width),
                                               y = as.double(1:map.width))),
                         to = c(0, 1))
# Distance matrix
D = as.matrix(dist(coords))

## exponential covariance function
exp_covar <- function(dst, rng){
  covar = exp(-dst/rng) # exponential covariance function
  return(covar)
}

# Calculate covariance matrix
Sigma = diag(rep(nug, npix)) + ((1 - nug) * exp_covar(dst = D, rng = r))


del = t(mvtnorm::rmvnorm(ntime, sigma = Sigma))


# temporal errors
eps = t(apply(X = del, MARGIN = 1, function(e){
  arima.sim(n = ntime, model = list(ar = AR, ma = 0, innov = e))
}))


data.frame(coords, del) %>%
  reshape2::melt(id.vars = c("x", "y")) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_tile() +
  facet_grid(variable~.)

data.frame(coords, eps) %>%
  reshape2::melt(id.vars = c("x", "y")) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_tile() +
  facet_grid(variable ~ .)

# spatial errors
delta = scale(Sigma %*% t(matrix(rnorm(npix*ntime, sd = .05), ncol = npix, nrow = ntime)))

delta = (Sigma %*% t(replicate(npix, rnorm(ntime))))

epsillon = t(apply(X = delta, MARGIN = 1, function(e){
  arima.sim(n = ntime, model = list(ar = AR, ma = 0, innov = e))
}))

data.frame(coords, delta) %>%
  reshape2::melt(id.vars = c("x", "y")) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_tile() +
  facet_grid(variable~.)

data.frame(coords, epsillon) %>%
  reshape2::melt(id.vars = c("x", "y")) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_tile() +
  facet_grid(variable ~ .)
