context("test-basic.R -- Basic Blip CDF Test")

# DGD function
gendata.blip=function(n, d, g0, Q0){
  # d=1
  # n=10
  var_names = paste0("W", 1:d)
  W_list = lapply(1:d, FUN = function(i) rnorm(n))
  names(W_list) = var_names
  pscores = do.call(g0, W_list)
  A=rbinom(n,1,pscores)
  AW_list = W_list
  AW_list$A = A
  Q = do.call(Q0, AW_list)
  Y=rbinom(n,1,Q)
  AW_list1 = AW_list0 = AW_list
  AW_list1$A = rep(1,n)
  AW_list0$A = rep(0,n)
  blip = do.call(Q0, AW_list1) - do.call(Q0, AW_list0)
  df = cbind(A = A, as.data.frame(W_list), Y = Y)
  return(list(df=df,blip=blip))
}


# Define th DGP functions for SCM
g0 = function(W1) plogis(.2 + .2*W1)
Q0 = function(A,W1) plogis(A + 2.5*A*W1 + W1)

Q_mean <- function(task) {
  W1 <- task$data$W1
  A <- task$data$A
  result <- Q0(A,W1)
  return(result)
}

g_dens <- function(task){
  W1 <- task$data$W1
  A <- task$data$A
  gA1 <- g0(W1)
  result <- ifelse(A==1, gA1, 1-gA1)
  return(result)
}

# make a grid of points at which to evaluate the CDF
m = -0.195
M = 0.319
cdf_points = seq(m, M, .01)

# define kernel

kernel = make_kernel(degree = 4, R = 2)
n = 1000
bw = n^-.2

# generate tmle task
data <- gendata.blip(1000, d = 1, g0, Q0)

node_list <- list(W="W1", A="A", Y="Y")

# make tmle_task
npsem <- list(
  define_node("W", node_list$W),
  define_node("A", node_list$A, c("W")),
  define_node("Y", node_list$Y, c("A", "W"))
)


tmle_task <- make_tmle3_Task(data$df, npsem)

# build likelihood object
# use known likelihoods
factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_known, "A", density_fun = g_dens),
  define_lf(LF_known, "Y", mean_fun = Q_mean, type = "mean")
)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)
likelihood$get_likelihoods(tmle_task)

# generate targeted likelihood
updater <- tmle3_Update$new(maxit=1e5, one_dimensional = TRUE, constrain_step = FALSE, delta_epsilon = 1e-2, verbose = TRUE)
targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)

# define parameter
intervention_treatment <- define_lf(LF_static, "A", value = 1)
intervention_control <- define_lf(LF_static, "A", value = 0)
bCDF <- define_param(Param_blipCDF, targeted_likelihood, intervention_treatment, intervention_control, kernel, cdf_points, bw)

# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, list(bCDF), updater)
estimates <- tmle_fit$estimates

# extract results
cdf_summary <- summary_from_estimates(tmle_task, estimates, simultaneous_ci = TRUE)
cdf_summary$blip <- cdf_points

library(ggplot2)
ggplot(cdf_summary, aes(x=blip, y= psi_transformed, ymin=lower_transformed,ymax=upper_transformed))+
  geom_line()+theme_bw()+ylab("P(B<t)")+xlab("t")+geom_ribbon(alpha=0.5)
