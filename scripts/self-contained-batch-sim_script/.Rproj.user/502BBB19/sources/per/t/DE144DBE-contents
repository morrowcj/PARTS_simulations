
# X(t) = lambda*land + S1*r1 + (alpha + beta*land + S2*r2 )*t + (gamma + S3*r3)*temp(t) e(t)

# spherical tapering function
taper_covar <- function(dst, rng){
  ts = ifelse(dst > rng, 0, ((1 - dst/rng)^2) * (1 + dst/(2 * rng)) )
  return(ts)
}
require(dplyr)
require(scales)
require(parallel)
require(geosphere)
require(remotePARTS)

## Simulation Function
sim_PARTs <- function(
  ## dimensions
  map.width = 48, # number of columns (in pixels)
  ntime = 20, # number of time points
  ## Model params
  lambda = 0.1, # base effect of land
  alpha = .05,  # base effect of time
  beta = .05, # effect of land:time
  gamma = .05, # base effect of temperature
  ## Variation params
  S1 = 0, # Intercept variation effect
  S2 = 0, # Time variation effect
  S3 = 0, # Temp variation effect
  r1_f = "r1_gen",
  r2_f = "r1_gen",
  r3_f = "r1_gen",
  ## spatial params
  r = .01, # range
  fitr.n = 1000, # number of samples to estimate R with
  nug = 0, # nugget
  AR = 0.2, # AR parameter for e(t)
  MA = 0, # MA parameter for e(t)

  AR_W = .2, # AR parameter for e_w(t)
  MA_W = 0, # MA parameter for e_w(t)

  land.grid = 4, # subdivision of land
  ## Climate paramters
  temp_f = "sin_temp", # function to generate base climate data from x, y coords
  temp.r = 1, # spatial range of temperature autocorrelation
  ## simulation params
  par.sims = 1,
  covar_f = "taper_covar", # covariance function
  debug = TRUE,
  verbose = TRUE,
  ## fixed variables among simulations
  fixed_temp = FALSE, #fixed climate variable?
  fixed_eps = FALSE, #fixed error term?
  fixed_r1 = FALSE, #fixed r1?
  fixed_r2 = FALSE, #fixed r2?
  fixed_r3 = FALSE, #fixed r3?

  t3_errors = FALSE, # should errors be non-guassian?
  err_sd = 1, # effect of the error term as in: e ~ N(0, err_sd*Sigma)

  npart = 10,
  partsize = NA,
  ncores = 8,

  GLS.resp = "time", # which coefficient to test? c("(Intercept)", "time", "z")
  int.only.form = FALSE # intercept only GLS model?
){


  if(debug){cat("debug mode\n")}

  ## Error Handling ----
  # check map width
  if (any(map.width %% 2 != 0) | any(map.width %% 4 != 0) | any(map.width %% 8 !=0)) {
    stop("map.width must be divisible by 2, 4, and 8")
  }

  call = match.call()

  ## Setup ----
  # sim parameters
  sim.table = expand.grid(map.width = map.width, lambda = lambda, alpha = alpha,
                          beta = beta, gamma = gamma, r = r, nug = nug,
                          AR = AR, MA = MA, AR_W = AR_W, MA_W = MA_W,
                          land.grid = land.grid, S1 = S1, S2 = S2, S3 = S3,
                          temp.r = temp.r, sim.no = seq_len(par.sims)) %>%
    arrange(map.width, land.grid, nug, r, temp.r)

  # number of sims
  nsims = nrow(sim.table)

  sim.table$simID <- seq_len(nsims)

  if(debug & nsims >= 3){nsims = 3} ## only do 3 loops max in debug mode

  # set starting values to character string
  map.width = land.grid = nug = r =  temp.r = "start"

  # match functions
  func = match.fun(covar_f) # covariance function
  temp_func = match.fun(temp_f) # climate generation function
  r1_func = match.fun(r1_f)
  r2_func = match.fun(r2_f)
  r3_func = match.fun(r3_f)

  ## Setup output ----
  # all_X = matrix(NA, nrow = npix*nsims, ncol = ntime)
  # all_land = matrix(NA, nrow = npix*nsims, ncol = 1)
  # all_Temp = matrix(NA, nrow = npix*nsims, ncol = ntime)
  GLS_outlist = list()

  ##
  # Sim Loop
  ##
  pb = txtProgressBar(min = 0, max = nsims, style = 3)

  for(sim.i in seq_len(nsims)){
    if (debug | verbose) {
      cat("sim", sim.i,"\n")
    }

    if (!debug & !verbose) {
        setTxtProgressBar(pb, sim.i)
      }

    ## Setup Map ----
    if (map.width != sim.table[sim.i, "map.width"]) { # build map things if they've changed
      map.width = sim.table[sim.i, "map.width"]
      size.changed = TRUE # indicate that the map has changed this loop
      # pixels in the map
      npix = map.width ^ 2
      # coordinates of each pixel
      coords = scales::rescale(as.matrix(expand.grid(x = as.double(1:map.width),
                                                     y = as.double(1:map.width))),
                               to = c(0, 1))

      # ## alternative way: (09-Sept-2021)
      # coords = cbind(x = rep(seq(0, 1, length.out = map.width), times = map.width),
      #                y = rep(seq(0, 1, length.out = map.width), each = map.width))

      # Distance matrix
      D = as.matrix(dist(coords))
    } else{
      size.changed = FALSE
    }

    ## Calculate covariance matrix ----
    if (nug != sim.table[sim.i, "nug"] | r != sim.table[sim.i, "r"] | size.changed){ # calculate Sigma when needed
      # set paramters
      nug = sim.table[sim.i, "nug"]; r = sim.table[sim.i, "r"]
      # covar mat
      Sigma = diag(rep(nug, npix)) + ((1 - nug) * func(dst = D, rng = r))
    }

    ## Establish land classes ----
    if (land.grid != sim.table[sim.i, "land.grid"] | size.changed) { # build land grid if they've changed
      land.grid = sim.table[sim.i, "land.grid"]
      if (land.grid == 1) {# all one land class
        pattern = matrix(1, ncol = 8, nrow = 8)
      } else if (land.grid == 2) { # 2 classes each 1/2 width of map
        pattern <- matrix(c(rep(c(1,0), each = 4, times = 4),
                            rep(c(0, 1), each = 4, times = 4)),
                          ncol = 8, byrow = TRUE)
      } else if (land.grid == 4) { # 2 classes each 1/4 width of map
        pattern <- matrix(rep(c(rep(c(1,0), each = 2, times = 4),
                                rep(c(0,1), each = 2, times = 4)), times = 2),
                          ncol = 8, byrow = TRUE)
      } else if (land.grid == 8) { # 2 classes 1/8 width of map
        pattern <-  matrix(rep(c(rep(c(1,0), times = 4),
                                 rep(c(0,1), times = 4)), times = 4),
                           ncol = 8, byrow = TRUE)
      } else { # error handling
        stop("land.grid must be one of 1, 2, 4, or 8")
      }
      # land pattern
      land = kronecker(pattern, matrix(1, map.width/8, map.width/8))
    }

    ## Setup climate Variable ----
    if (temp.r != sim.table[sim.i, "temp.r"] | size.changed){ # calculate climate covariance when neccessary
      temp.r = sim.table[sim.i, "temp.r"]
      # baseline temperature pattern
      basetemp = temp_func(coords[, "x"], coords[, "y"], npix)
      # get parameters
      temp.r = sim.table[sim.i, "temp.r"]
      # spatial covariance
      temp.Sig = diag(rep(0, npix)) + ((1-0) * func(dst = D, rng = temp.r))
    }
    AR_W = sim.table[sim.i, "AR_W"]; MA = sim.table[sim.i, "MA_W"]
    if (fixed_temp) {
      # Fixed temperature
      if (size.changed) { # calculate temp only once for a given map size
        # spatial errors
        if (ntime == 1){
          temp.del = scale(temp.Sig %*% matrix(rnorm(ntime), nrow = npix))
        } else {
          temp.del = scale(temp.Sig %*% t(replicate(npix, rnorm(ntime))))
        }
        # temporal errors
        if(ntime == 1){
          temp.eps = matrix(apply(X = temp.del, MARGIN = 1, function(e){
            arima.sim(n = ntime, model = list(ar = AR_W, ma = MA_W, innov = e))
          }), ncol = ntime)
        } else {
          temp.eps = t(apply(X = temp.del, MARGIN = 1, function(e){
            arima.sim(n = ntime, model = list(ar = AR_W, ma = MA_W, innov = e))
          }))
        }
        # spatio-temporal climate variable
        Temp = replicate(ntime, as.vector(basetemp)) +
          matrix(rep(seq_len(ntime)*.1, each = npix), ncol = ntime) +
          scale(temp.eps)
        colnames(Temp) <- paste("Temp", 1:ntime, sep = ".")
      }
    } else { # calculate temp for each loop
      # Variable temperature
      # spatial errors
      if (ntime == 1){
        temp.del = scale(temp.Sig %*% matrix(rnorm(ntime), nrow = npix) )
      } else {
        temp.del = scale(temp.Sig %*% t(replicate(npix, rnorm(ntime))) )
      }
      # temporal errors
      if(ntime == 1){
        temp.eps = matrix(apply(X = temp.del, MARGIN = 1, function(e){
          arima.sim(n = ntime, model = list(ar = AR_W, ma = MA_W, innov = e))
        }), ncol = ntime)
      } else {
        temp.eps = t(apply(X = temp.del, MARGIN = 1, function(e){
          arima.sim(n = ntime, model = list(ar = AR_W, ma = MA_W, innov = e))
        }))
      }
      # spatio-temporal climate variable
      Temp = replicate(ntime, as.vector(basetemp)) +
        matrix(rep(seq_len(ntime)*.1, each = npix), ncol = ntime) +
        scale(temp.eps)
      colnames(Temp) <- paste("Temp", 1:ntime, sep = ".")
    }

    ## Calculate spatio-temporal errors ----
    AR = sim.table[sim.i, "AR"]; MA = sim.table[sim.i, "MA"]
    if (fixed_eps) {
      if (size.changed){ # calculate epsilon only once for each map size
        # spatial errors
        if (ntime == 1){
          if(t3_errors){
            delta = scale(Sigma %*% matrix(rt(ntime, 3), nrow = npix))
          } else {
            delta = scale(Sigma %*% matrix(rnorm(ntime), nrow = npix))
          }
        } else {
          if(t3_errors){
            delta = scale(Sigma %*% matrix(rt(ntime, 3), nrow = npix))
          } else {
            delta = scale(Sigma %*% t(replicate(npix, rnorm(ntime))))
          }
        }
        # temporal errors
        if(ntime == 1){
          eps = matrix(apply(X = delta, MARGIN = 1, function(e){
            arima.sim(n = ntime, model = list(ar = AR, ma = MA, innov = e))
          }), ncol = ntime)
        } else {
          eps = t(apply(X = delta, MARGIN = 1, function(e){
            arima.sim(n = ntime, model = list(ar = AR, ma = MA, innov = e))
          }))
        }
      }
    } else { # calculate epsilon for each iteration
      # spatial errors
      if (ntime == 1){
        delta = scale(Sigma %*% matrix(rnorm(ntime), nrow = npix))
      } else {
        delta = scale(Sigma %*% t(replicate(npix, rnorm(ntime))))
      }
      # temporal errors
      if(ntime == 1){
        eps = matrix(apply(X = delta, MARGIN = 1, function(e){
          arima.sim(n = ntime, model = list(ar = AR, ma = MA, innov = e))
        }), ncol = ntime)
      } else {
        eps = t(apply(X = delta, MARGIN = 1, function(e){
          arima.sim(n = ntime, model = list(ar = AR, ma = MA, innov = e))
        }))
      }
    }

    ## Setup umeasured random variables ----
    # Intercept
    S1 = sim.table[sim.i, "S1"]
    if (fixed_r1) {
      if (size.changed) {
        r1 = r1_func(coords[, "x"], coords[, "y"], npix)
        Sr1 = replicate(ntime, as.vector(S1 * r1))
      }
    } else {
      r1 = r1_func(coords[, "x"], coords[, "y"], npix)
      Sr1 = replicate(ntime, as.vector(S1 * r1))
    }
    # Time
    S2 = sim.table[sim.i, "S2"]
    if (fixed_r2) {
      if (size.changed) {
        r2 = r2_func(coords[, "x"], coords[, "y"], npix)
        Sr2 = replicate(ntime, as.vector(S2 * r2))
      }
    } else {
      r2 = r2_func(coords[, "x"], coords[, "y"], npix)
      Sr2 = replicate(ntime, as.vector(S2 * r2))
    }
    # Climate
    S3 = sim.table[sim.i, "S3"]
    if (fixed_r3) {
      if (size.changed) {
        r3 = r3_func(coords[, "x"], coords[, "y"], npix)
        Sr3 = replicate(ntime, as.vector(S3 * r3))
      }
    } else {
      r3 = r3_func(coords[, "x"], coords[, "y"], npix)
      Sr3 = replicate(ntime, as.vector(S3 * r3))
    }

    ## Build X ----
    # Intercept effect
    lambda = sim.table[sim.i, "lambda"]
    lamb.land = matrix(rep(lambda * land, times = ntime), nrow = npix, ncol = ntime)
    int.eff =  lamb.land + Sr1
    # Time effect
    alpha = sim.table[sim.i, "alpha"]
    beta = sim.table[sim.i, "beta"]
    if(ntime == 1){
      alpha.t = matrix(alpha, nrow = npix, ncol = ntime)  * ntime
      beta.land.t = matrix(beta, nrow = npix, ncol = ntime) * ntime
      Sr2.t = Sr2 * ntime
    } else {
      alpha.t = replicate(ntime, rep(alpha, npix)) * t(replicate(npix, seq_len(ntime)))
      beta.land.t = replicate(ntime, rep(beta, npix)*as.vector(land)) * t(replicate(npix, seq_len(ntime)))
      Sr2.t = Sr2 * t(replicate(npix, seq_len(ntime)))
    }
    time.eff = alpha.t + beta.land.t + Sr2.t
    # Climate effect
    gamma = sim.table[sim.i, "gamma"]
    temp.eff = gamma*Temp + Sr3
    # X
    X = int.eff + time.eff + temp.eff + err_sd*eps
    colnames(X) <- paste("X", 1:ntime, sep = ".")

    full.dat = data.frame(x.coord = coords[,"x"], y.coord = coords[,"y"],
                          land = factor(as.vector(land)), X, Temp)


    require(remotePARTS)

    if(ntime == 1){
      full.dat$C.coef = X
    } else if (ntime == 2) {
      stop("ntime cannot be 2. AR needs at least 3 time points")
    } else {

      # response = GLS.resp
      # fm.AR = fitAR.map(X, t = seq_len(ntime), Z = Temp, ret_int.coef = TRUE, ret_AR.par = TRUE)
      # full.dat$C.coef = fm.AR[[response]][,"Est"] # C = alpha + beta + S2

      if (debug | verbose) {cat("AR step\n")}
      clst <- parallel::makeCluster(ncores) # initiate cluster: 8 cores
      # clusterExport(clst, varlist = c("ntime", "GLS.resp"))
      full.dat$C.coef <- parApply(cl = clst, X = cbind(X, Temp), MARGIN = 1,
                                  FUN = function(x, nt = ntime){
                                    library(remotePARTS)
                                    y = x[1:nt]
                                    z = x[(nt+1):(2*nt)]
                                    time = seq_len(nt)
                                    fitAR(y ~ time + z)$coef[GLS.resp, "Est"]
                                  }
      )
      stopCluster(clst)

    }

    ## estimate r ----
    # r.hat <- tryCatch(fit_spatialcor(X, t = seq_len(ntime), method = "exponential",
    #                                  fit.n = 1000, location = coords,
    #                                  r.start = .05)$spatialcor,
    #                   error = function(e){
    #                     ## What do do if there's an error
    #                     cat("error estimating r\n")
    #                     return(.1) # for now, just return a useable r
    #                   })

    if (ntime == 1){
      r.hat = r/2.5
      if (debug | verbose) {cat("assuming r is known to be",r.hat,"\n")}
    } else {
    if (debug | verbose) {cat("estimating r\n")}
    samp.int = sample.int(npix, fitr.n, replace = FALSE)
    clst <- parallel::makeCluster(ncores) # initiate cluster: 8 cores
    # clusterExport(clst, varlist = c("ntime", "GLS.resp"))
    X.resids = parApply(cl = clst, X = cbind(X[samp.int, ], Temp[samp.int, ]), MARGIN = 1,
                        FUN = function(x, nt = ntime){
                          library(remotePARTS)
                          y = x[1:nt]
                          z = x[(nt+1):(2*nt)]
                          time = seq_len(nt)
                          fitAR(y ~ time + z)$resids
                        }
    )
    stopCluster(clst)

    cor.resid = cor(X.resids) # as a matrix
    vec_cor = cor.resid[upper.tri(cor.resid, diag = TRUE)] # as a vector
    Dsub = D[samp.int, samp.int]
    # turn the distance matrix into a vector and scale it to between 0, and 1 for simpicity
    vec_dist = Dsub[upper.tri(Dsub, diag = TRUE)] # as a vector
    W = data.frame(dist = vec_dist + .0001, cor = vec_cor) # combine into data frame
    # try to estimate r
    # fit = nls(cor ~ exp(-(dist/r)), data = W, start = list(r = 1), algorithm = "port")
    # r.hat = coef(fit)[1]

    r.hat <- tryCatch(expr = {
      fit = nls(cor ~ exp(-(dist/r)), data = W, start = list(r = .1), algorithm = "port")
      coef(fit)[1]
    },
    error = function(e){
      print(e)
      cat("error estimating r. Using r = 0.004\n")
      return(.01/2.5)
    })

    }


    ## GLS Step ----

    if (debug | verbose) {cat("GLS\n")}

    part_M = sample_partitions(npix, partsize = partsize, npart = npart, verbose = FALSE)
    part_dims = dim(part_M)

    if(int.only.form){
      form = "C.coef ~ 1"
    } else {
      form = "C.coef ~ 1 + land"
    }

    GLS = fitGLS.partition.mc(part_f = "part_data", dist_f = "dist_tmp",
                              V.meth = "exponential", spatcor = r.hat,
                              partsize = part_dims[1], npart = part_dims[2],
                              mincross = 6, ncores = ncores,
                              debug = debug,
                              part_form = form,
                              part_df = full.dat,
                              part_mat = part_M,
                              part_locvars = c("x.coord", "y.coord"))

    t.pval = cor_t.test(GLS)
    chisq.pval = cor_chisq.test(GLS)

    GLS_outlist[[sim.i]] <- list(part_coefs = GLS$part.stats$coefficients,
                                 part_modstats = GLS$part.stats$mod.stats,
                                 cross_rSS = GLS$cross.stats$cross.SS,
                                 cross_coef = GLS$cross.stats$rcoefs,
                                 overall_stats = GLS$overall.stats,
                                 t.pval = t.pval,
                                 chisq.pval = chisq.pval)

    other.info = c(GLS.resp = GLS.resp, int.only = int.only.form, par.sims = par.sims, nsims = nsims)

  }
  if (debug) {
    return(list(mod.info = list(call = call, other = other.info),
                GLS = GLS_outlist,
                sim.table = sim.table,
                data = full.dat,
                delta = delta,
                eps = eps,
                r1 = r1,
                # AR = fm.AR,
                r.hat = r.hat
    ))
  } else {
    return(list(mod.info = list(call = call, other = other.info),
                int.only.mod = int.only.form,
                GLS = GLS_outlist,
                sim.table = sim.table))
  }
}

# distance function
dist_tmp <- function(coords, coords2 = NULL){
  if(is.null(coords2)){
    return(geosphere::distm(coords)/100000)
  } else {
    return(geosphere::distm(coords, coords2)/100000)
  }
}

# exponential covariance function
exp_covar <- function(dst, rng){
  covar = exp(-dst/rng) # exponential covariance function
  return(covar)
}

# spherical tapering function
taper_covar <- function(dst, rng){
  ts = ifelse(dst > rng, 0, ((1 - dst/rng)^2) * (1 + dst/(2 * rng)) )
  return(ts)
}

# 1D sin wave
sw = function(w, A = 1, T = 50/3, phi = 0){
  A*sin(((2*pi) / T)*w + phi)
}

# 2D sin wave
sin_2d = function(x, y, A = 1, Tx = .6, Ty = .6, px = 45, py = 45){
  A * sin( (2*pi/Tx)*x + px + (2*pi/Ty)*y + py )
}

# Base temperature from 2D sin wave
sin_temp <- function(coord.x, coord.y, npix){
  scale(sin_2d(coord.x, coord.y, Tx = -5, Ty = -3, px = 45, py = 45) + rnorm(npix, sd = 1/4))
}

# function to generate random variables
r1_gen <- function(coord.x, coord.y, npix){
  (scale(sin_2d(coord.x, coord.y, npix, Tx = .6, Ty = .6, px = 45, py = 45)) + rnorm(npix))
}
