#!/usr/bin/env Rscript

# load libraries
require(reshape2)
library(getopt)

##
# Parse Arguments ----
##

## Parse command arguments #0: no args, 1: mandatory args, 2: optional args
spec = matrix(c(
  'help'    , 'h', 0, "logical",
  'debug'   , 'u', 0, "logical",
  "cores"   , "c", 1, "integer",
  "testrun" , "t", 0, "logical",
  "notrun"  , "n", 0, "logical",
  "notoneA" , "i", 0, "logical",
  "notoneB" , "j", 0, "logical",
  "notoneC" , "k", 0, "logical",
  "nottwo"  , "l", 0, "logical",
  "notthree", "f", 0, "logical",
  "verbose" , "v", 0, "logical",
  "savesims", "s", 0, "logical"),
  ncol = 4, byrow = TRUE)
opt = getopt(spec)

## help was asked for - immediately exit
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage = TRUE))
  q(status=1);
}

## defaulf values
if (is.null(opt$debug)) {opt$debug = FALSE}
if (is.null(opt$verbose)) {opt$verbose = FALSE}
if (is.null(opt$testrun)) {opt$testrun = FALSE}
if (is.null(opt$cores)) {opt$cores = 4}
if (is.null(opt$notrun)) {opt$notrun = FALSE}
if (is.null(opt$notoneA)) {opt$notoneA = FALSE}
if (is.null(opt$notoneB)) {opt$notoneB = FALSE}
if (is.null(opt$notoneC)) {opt$notoneC = FALSE}
if (is.null(opt$nottwo)) {opt$nottwo = FALSE}
if (is.null(opt$notthree)) {opt$notthree = FALSE}
if (is.null(opt$savesims)) {opt$savesims = FALSE}



## debug mode
if (opt$debug) {
  cat("debug mode\n")
}

##
# Setup ----
##

sims.default = 200 # change this number to do more sims.

## Function to get root dir of this script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

# save the root directory
if (rstudioapi::isAvailable()) { # within Rstudio
  base.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
} else { # within Rterm or Rscript
  base.dir = dirname(thisFile())
}

# source the simulation function file
source(file.path(base.dir, "simulation-functions_morrowcj.R"))



##
# Study 2-1a: Effect of random spatial autocorrelation ----
##
# # Full Model:
# X(t) = lambda*land + S1*r1 + (alpha + beta*land + S2*r2 )*t + (gamma + S3*r3)*temp(t) e(t)
#
# # Temporal error:
# e(t) ~ rho*e(t - 1) + d(t)
#
# # Spatial error:
# d ~ N(0, Sigma)
#
# # Spatial covariance:
# Sigma = f(D, theta, nug)
#
# # Sim Model:
# X(1) = lambda*land + e(1)
# - lamda = 0.2
# - S1 = 0
# - alpha = 0
# - beta = 0
# - S2 = 0
# - gamma = 0
# - S3 = 0
# - theta = {0, 0.05, 0.25}


if (!opt$notrun & !opt$notoneA) {

  cat("\nSim study 2-1a\n")

  # if --testrun is supplied, make data small
  if (opt$testrun) {map.width = c(32)} else {map.width = c(104)}
  if (opt$testrun) {sims = 4} else {sims = 500}
  if (opt$testrun) {partsize = 100} else {partsize = 2000}
  if (opt$testrun) {fitr.n = 100} else {fitr.n = 1000}


  if(!dir.exists(file.path(base.dir,"sim-data"))){
    dir.create(file.path(base.dir, "sim-data"))
  }

  sim2.1a_time <- system.time(expr = {
    sim2.1a = sim_PARTs(map.width = map.width,
                     ntime = 1, # only spatial
                     land.grid = 4, # granularity 4
                     ## model params
                     lambda = 0.2, # fixed effect of L2
                     S1 = 0, fixed_r1 = TRUE, # no fixed spatial autocor
                     alpha = 0, beta = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                     gamma = 0, S3 = 0, fixed_r3 = TRUE, fixed_temp = TRUE, # no climate variable
                     r = c(1e-9, 0.05, 0.25), # variable range
                     nug = 0, # no nugget
                     ## parameter of interest
                     GLS.resp = "(Intercept)", # test lambda
                     ## logistiscs pars
                     t3_errors = FALSE,
                     fitr.n = fitr.n, partsize = partsize, par.sims = sims, npart = NA,
                     verbose = opt$debug, debug = opt$debug, ncores = opt$cores
                     )
  })

  if(opt$debug){cat("sim 2-1a time:", sim2.1a_time[1], "\n")}
  if (!opt$testrun){
    save(sim2.1a, sim2.1a_time,
         file = file.path(base.dir,"sim-data/2-01a_PARTS_sim-study.RData"))
  }
} else{
  load(file.path(base.dir, "sim-data/2-01a_PARTS_sim-study.RData"))
}

##
# Study 2-1b.1: non-guassian errors ----
##
#
# # Full Model:
# X(t) = lambda*land + S1*r1 + (alpha + beta*land + S2*r2 )*t + (gamma + S3*r3)*temp(t) e(t)
#
# # Temporal error:
# e(t) ~ rho*e(t - 1) + d(t)
#
# # Spatial error:
# d ~ Sigma*T(3)
#
# # Spatial covariance:
# Sigma = f(D, theta, nug)
#
# # Sim Model:
# X(1) = lambda*land + e(1)
# - lamda = {0, 0.2}
# - S1 = 0
# - alpha = 0
# - beta = 0
# - S2 = 0
# - gamma = 0
# - S3 = 0
# - theta = {0, 0.05, 0.25}

if (!opt$notrun & !opt$notoneB) {

  cat("\nSim study 2-1b\n")

  # if --testrun is supplied, make simulation data small
  if (opt$testrun) {map.width = c(32)} else {map.width = c(104)}
  if (opt$testrun) {
    sims.1 = 4
    sims.2 = 4
  } else {
    sims.1 = 500
    sims.2 = 1000
    }
  if (opt$testrun) {partsize = 100} else {partsize = 2000}
  if (opt$testrun) {fitr.n = 100} else {fitr.n = 1000}

  sim2.1b_time <- system.time(expr = {
    # 1b1: 500 sims for variable r
    cat("   i\n")
    sim2.1b1 = sim_PARTs(map.width = map.width,
                        ntime = 1, # only spatial
                        land.grid = 4, # granularity 4
                        ## model params
                        lambda = 0.2, # fixed effect of L2
                        S1 = 0, fixed_r1 = TRUE, # no fixed spatial autocor
                        alpha = 0, beta = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                        gamma = 0, S3 = 0, fixed_r3 = TRUE, fixed_temp = TRUE, # no climate variable
                        r = c(1e-9, 0.05, 0.25), # variable range
                        nug = 0, # no nugget
                        ## parameter of interest
                        GLS.resp = "(Intercept)", # test lambda
                        ## logistiscs pars
                        t3_errors = TRUE,
                        fitr.n = fitr.n, partsize = partsize, par.sims = sims.1, npart = NA,
                        verbose = opt$debug, debug = opt$debug, ncores = opt$cores,
    )
    # 1b2: 1000 sims for type I errors
    cat("\n   ii\n")
    sim2.1b2 = sim_PARTs(map.width = map.width,
                         ntime = 1, # only spatial
                         land.grid = 4, # granularity 4
                         ## model params
                         lambda = 0, # fixed effect of L2
                         S1 = 0, fixed_r1 = TRUE, # no fixed spatial autocor
                         alpha = 0, beta = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                         gamma = 0, S3 = 0, fixed_r3 = TRUE, fixed_temp = TRUE, # no climate variable
                         r = c(0.05), # fixed range
                         nug = 0, # no nugget
                         ## parameter of interest
                         GLS.resp = "(Intercept)", # test lambda
                         ## logistiscs pars
                         t3_errors = TRUE,
                         fitr.n = fitr.n, partsize = partsize, par.sims = sims.2, npart = NA,
                         verbose = opt$debug, debug = opt$debug, ncores = opt$cores,
    )
  })

  if(opt$debug){cat("sim 2-1b time:", sim2.1b_time[1], "\n")}
  if (!opt$testrun | opt$savesims) {
    save(sim2.1b1, sim2.1b2, sim2.1b_time,
         file = file.path(base.dir, "sim-data/2-01b_PARTS_sim-study.RData"))
  }
} else{
  load(file.path(base.dir, "sim-data/2-01b_PARTS_sim-study.RData"))
}

##
# Study 2-1c ----
##
#
# # Full Model:
# X(t) = lambda*land + S1*r1 + (alpha + beta*land + S2*r2 )*t + (gamma + S3*r3)*temp(t) e(t)
#
# # Temporal error:
# e(t) ~ rho*e(t - 1) + d(t)
#
# # Spatial error:
# d ~ N(0, Sigma)
#
# # Spatial covariance:
# Sigma = f(D, theta, nug)
#
# # Sim Model:
# X(1) = lambda*land + e(1)
# - lamda = 0.2
# - S1 = sd(e);
# - R1 ~ 2d_sin(coords, npeaks); npeaks = {1^2, 2^2, 3^2}
# - alpha = 0
# - beta = 0
# - S2 = 0
# - gamma = 0
# - S3 = 0
# - theta = 0

if (!opt$notrun & !opt$notoneC) {
  cat("\nSim study 2-1c\n")

  # if --testrun is supplied, make simulation data small
  if (opt$testrun) {map.width = c(32)} else {map.width = c(104)}
  if (opt$testrun) {sims = 4} else {sims = 200}
  if (opt$testrun) {partsize = 100} else {partsize = 2000}
  if (opt$testrun) {fitr.n = 100} else {fitr.n = 1000}



  sim2.1c_time <- system.time(expr = {
    S1 = 0.5
    err_sd = 0.5
    # sin function with 1 peak
    one_peak <- function(coord.x, coord.y, npix, npeaks = 1^2){
      tx = ty = 1/(npeaks/2)
      -scale(sin_2d(coord.x, coord.y, Tx = tx, Ty = ty, px = 45, py = 45, A = npix))
    }
    # 1c1 - npeaks = 1^2
    cat("   i\n")
    sim2.1c1 = sim_PARTs(map.width = map.width,
                        ntime = 1, # only spatial
                        land.grid = 4, # granularity 4
                        ## model params
                        lambda = 0.2, # fixed effect of L2
                        S1 = S1, fixed_r1 = TRUE, # no fixed spatial autocor
                        alpha = 0, beta = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                        gamma = 0, S3 = 0, fixed_r3 = TRUE, fixed_temp = TRUE, # no climate variable
                        r = 1e-9, # variable range
                        nug = 0, # no nugget
                        err_sd = err_sd,
                        ## parameter of interest
                        GLS.resp = "(Intercept)", # test lambda
                        ## logistiscs pars
                        t3_errors = FALSE,
                        fitr.n = fitr.n, partsize = partsize, par.sims = sims, npart = NA,
                        verbose = opt$debug, debug = opt$debug, ncores = opt$cores,
                        r1_f = "one_peak"
    )
    # sin function with 4 peaks
    four_peaks <- function(coord.x, coord.y, npix, npeaks = 2^2){
      tx = ty = 1/(npeaks/2)
      -scale(sin_2d(coord.x, coord.y, Tx = tx, Ty = ty, px = 45, py = 45, A = npix))
    }
    # 1c2 - npeaks = 2^2
    cat("\n   ii\n")
    sim2.1c2 = sim_PARTs(map.width = map.width,
                        ntime = 1, # only spatial
                        land.grid = 4, # granularity 4
                        ## model params
                        lambda = 0.2, # fixed effect of L2
                        S1 = S1, fixed_r1 = TRUE, # no fixed spatial autocor
                        alpha = 0, beta = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                        gamma = 0, S3 = 0, fixed_r3 = TRUE, fixed_temp = TRUE, # no climate variable
                        r = 1e-9, # variable range
                        nug = 0, # no nugget
                        err_sd = err_sd,
                        ## parameter of interest
                        GLS.resp = "(Intercept)", # test lambda
                        ## logistiscs pars
                        t3_errors = FALSE,
                        fitr.n = fitr.n, partsize = partsize, par.sims = sims, npart = NA,
                        verbose = opt$debug, debug = opt$debug, ncores = opt$cores,
                        r1_f = "four_peaks"
    )
    # sin function with 6 peaks
    six_peaks <- function(coord.x, coord.y, npix, npeaks = 3^2){
      tx = ty = 1/(npeaks/2)
      -scale(sin_2d(coord.x, coord.y, Tx = tx, Ty = ty, px = 45, py = 45, A = npix))
    }
    # 1c3 - npeaks = 3^2
    cat("\n   iii\n")
    sim2.1c3 = sim_PARTs(map.width = map.width,
                        ntime = 1, # only spatial
                        land.grid = 4, # granularity 4
                        ## model params
                        lambda = 0.2, # fixed effect of L2
                        S1 = S1, fixed_r1 = TRUE, # no fixed spatial autocor
                        alpha = 0, beta = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                        gamma = 0, S3 = 0, fixed_r3 = TRUE, fixed_temp = TRUE, # no climate variable
                        r = 1e-9, # variable range
                        nug = 0, # no nugget
                        err_sd = err_sd,
                        ## parameter of interest
                        GLS.resp = "(Intercept)", # test lambda
                        ## logistiscs pars
                        t3_errors = FALSE,
                        fitr.n = fitr.n, partsize = partsize, par.sims = sims, npart = NA,
                        verbose = opt$debug, debug = opt$debug, ncores = opt$cores
    )
  })

  if(opt$debug){cat("sim 2-1c time:", sim2.1c_time[1], "\n")}
  if (!opt$testrun | opt$savesims) {
    save(sim2.1c1, sim2.1c2, sim2.1c3, sim2.1c_time,
         file = file.path(base.dir, "sim-data/2-01c_PARTS_sim-study.RData"))
  }
} else{
  load(file.path(base.dir, "sim-data/2-01c_PARTS_sim-study.RData"))
}

##
# Study 2a ----
##
#
# # Full Model:
# X(t) = lambda*land + S1*r1 + (alpha + beta*land + S2*r2 )*t + (gamma + S3*r3)*temp(t) e(t)
#
# # Temporal error:
# e(t) ~ rho*e(t - 1) + d(t)
#
# # Spatial error:
# d ~ N(0, Sigma)
#
# # Spatial covariance:
# Sigma = f(D, theta, nug)
#
# # Sim Model:
# X(t) = lambda*land + beta*land*t + e(t)
# - t = 1, ..., 30
# - lamda = 0.2
# - S1 = 0;
# - alpha = 0
# - beta = 1/30
# - S2 = 0
# - gamma = 0
# - S3 = 0
# - theta = {0, 0.05, 0.25}
# - nug = 0

if (!opt$notrun & !opt$nottwo) {
  cat("\nSim study 2-2a\n")

  # if --testrun is supplied, make simulation data small
  if (opt$testrun) {map.width = 32} else {map.width = 104}
  if (opt$testrun) {sims.a = 4} else {sims = 200}
  if (opt$testrun) {partsize = 100} else {partsize = 2000}
  if (opt$testrun) {fitr.n = 100} else {fitr.n = 1000}

  ## S3: time trend effects
  sim2.2_time <- system.time(expr = {
    sim2.2 = sim_PARTs(map.width = map.width,
                         ntime = 30, # 30 year series
                         land.grid = 4, # granularity 4
                         ## model params
                         lambda = 0.2, # fixed effect of L2
                         S1 = 0, fixed_r1 = TRUE, # no fixed spatial autocor
                         alpha = 0, S2 = 0, fixed_r2 = TRUE, # no time intercept or R
                         beta = 1/30,  # fixed time trend when L=1
                         gamma = 0, fixed_temp = TRUE, S3 = 0, fixed_r3 = TRUE, # no climate
                         r = c(0, 0.05, 0.25), # variable range
                         nug = 0, # no nugget
                         err_sd = 1,
                         AR = 0.4,
                         ## parameter of interest
                         GLS.resp = "time", # test beta
                         ## logistiscs pars
                         t3_errors = FALSE,
                         fitr.n = fitr.n, partsize = partsize, par.sims = sims, npart = NA,
                         verbose = opt$debug, debug = opt$debug, ncores = opt$cores
    )
  })

  if(opt$debug){cat("sim 2-2 time:", sim2.2_time[1], "\n")}
  if (!opt$testrun | opt$savesims){
    save(sim2.2, sim2.2_time,
         file = file.path(base.dir, "sim-data/2-02_PARTS_sim-study.RData"))
  }
} else{
  load(file.path(base.dir, "sim-data/2-02_PARTS_sim-study.RData"))
}

##
# Study 2 ----
##
#
# # Full Model:
# X(t) = lambda*land + S1*r1 + (alpha + beta*land + S2*r2 )*t + (gamma + S3*r3)*temp(t) e(t)
#
# # Temporal error:
# e(t) ~ rho*e(t - 1) + d(t)
#
# # Spatial error:
# d ~ N(0, Sigma)
#
# # Spatial covariance:
# Sigma = f(D, theta, nug)
#
# # Sim Model:
# X(t) = lambda*land + gamma*temp(t) + e(t)
# - t = 1, ..., 30
# - lamda = 0.2
# - S1 = 0;
# - alpha = 0
# - beta = 0
# - S2 = 0
# - gamma = 1
# - S3 = 0
# - theta = {0, 0.05}
# - theta_w = {0, 0.04}
# - nug = 0

if (!opt$notrun & !opt$notthree) {
  cat("\nSim study 2-3\n")

  # if --testrun is supplied, make simulation data small
  if (opt$testrun) {map.width = c(32)} else {map.width = c(104)}
  if (opt$testrun) {sims = 4} else {sims = 200}
  if (opt$testrun) {partsize = 100} else {partsize = 2000}
  if (opt$testrun) {fitr.n = 100} else {fitr.n = 1000}

  ## S3: climate effects
  sim2.3_time <- system.time(expr = {
    sim2.3 = sim_PARTs(map.width = map.width,
                        ntime = 30, # 30 year series
                        land.grid = 4, # granularity 4
                        ## model params
                        lambda = 0.2, # fixed effect of L2
                        S1 = 0, fixed_r1 = TRUE, # no fixed spatial autocor
                        beta = 0, alpha = 0, S2 = 0, fixed_r2 = TRUE, # no time trend
                        gamma = 1, fixed_temp = FALSE, # random climate varition
                        S3 = 0, fixed_r3 = TRUE, # no fixed climate variation
                        r = 0.25, # fixed range
                        nug = 0, # no nugget
                        err_sd = 1,
                        ## additional spatio-temporal params
                        AR = 0.4,
                        ## spatio-temporal climate params
                        AR_W = c(0, 0.4), # variable climate AR
                        temp.r = c(0, 0.25), # variable climate range
                        ## parameter of interest
                        GLS.resp = "z", # test beta
                        ## logistiscs pars
                        t3_errors = FALSE,
                        fitr.n = fitr.n, partsize = partsize, par.sims = sims, npart = NA,
                        verbose = opt$debug, debug = opt$debug, ncores = opt$cores
    )
  })


  if(opt$debug){cat("sim 2-3 time:", sim2.2_time[1], "\n")}
  if (!opt$testrun | opt$savesims){
    save(sim2.3, sim2.3_time,
         file = file.path(base.dir, "sim-data/2-03_PARTS_sim-study.RData"))
  }
} else{
  load(file.path(base.dir, "sim-data/2-03_PARTS_sim-study.RData"))

}

## End ----

if (! interactive()){
  quit(status = 0)
}
