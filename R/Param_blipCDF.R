#' Blip CDF
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
Param_blipCDF <- R6Class(
  classname = "Param_blipCDF",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, kernel, cdf_points, bandwidth, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      private$.kernel <- kernel
      private$.cdf_points <- cdf_points
      private$.bandwidth <- bandwidth
    },
    clever_covariates = function(tmle_task = NULL, cv_fold = -1) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      intervention_nodes <- names(self$intervention_list_treatment)
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes)
      
      # todo: think about collapse for multiple intervention nodes
      H_ATE <- (cf_pA_treatment - cf_pA_control) / pA

      cf_task_treatment <- self$cf_likelihood_treatment$cf_tasks[[1]]
      cf_task_control <- self$cf_likelihood_control$cf_tasks[[1]]

      EY1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, self$outcome_node)
      EY0 <- self$observed_likelihood$get_likelihoods(cf_task_control, self$outcome_node)
      
      B <- EY1 - EY0
      nn <- tmle_task$nrow
      

   
      
      HA <- vapply(self$cdf_points, FUN = function(x) {
        (1/self$bandwidth)*with(self$kernel, kern((B-x)/self$bandwidth, R=R, veck=veck))*H_ATE
      } ,FUN.VALUE= rep(1,nn))
    

      return(list(Y = unlist(HA, use.names = FALSE)))
    },
    estimates = function(tmle_task = NULL, cv_fold = -1) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      cf_task_treatment <- self$cf_likelihood_treatment$cf_tasks[[1]]
      cf_task_control <- self$cf_likelihood_control$cf_tasks[[1]]
      
      EY1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, self$outcome_node)
      EY0 <- self$observed_likelihood$get_likelihoods(cf_task_control, self$outcome_node)
      EYA <- self$observed_likelihood$get_likelihoods(tmle_task, self$outcome_node)
      
      B <- EY1 - EY0
      nn <- tmle_task$nrow
      
      
      
      int <- vapply(self$cdf_points, FUN = function(x0) {
        w = with(self$kernel, kern_cdf((B - x0)/self$bandwidth, R=R, veck=veck))
        return(w)
      } ,FUN.VALUE=rep(1,nn))
      psi <- apply(int, 2, mean)
      
      HA <- self$clever_covariates(tmle_task)$Y
      Y <- tmle_task$get_tmle_node("Y")
      # IC <- vapply(1:length(self$cdf_points),FUN = function(x) {HA[,x]*(Y-EYA)+int[,x]-psi[x]}, FUN.VALUE = rep(1,nn))
      IC <- HA*as.vector(Y-EYA)+int-rep(psi,each=nn)
      
      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("P(B < %0.3f)", self$cdf_points)
      return(param_form)
    },
    cf_likelihood_treatment = function() {
      return(private$.cf_likelihood_treatment)
    },
    cf_likelihood_control = function() {
      return(private$.cf_likelihood_control)
    },
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    kernel = function() {
      return(private$.kernel)
    },
    cdf_points = function() {
      return(private$.cdf_points)
    },
    bandwidth = function() {
      return(private$.bandwidth)
    },
    update_nodes = function() {
      return(self$outcome_node)
    }
  ),
  private = list(
    .type = "blipCDF",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .kernel = NULL,
    .cdf_points = NULL,
    .bandwidth = NULL
  )
)