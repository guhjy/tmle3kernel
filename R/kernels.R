
#' Generate a Polynomial Kernel
#' @param degree the degree of the polynomial for the kernel
#' @param R the domain of the kernel
#' @export 
make_kernel = function(degree, R){
  if (is.null(degree)) {
    kern = function(x, R, veck) 1/(2*R)*as.numeric(-R <= x & R >= x) 
    kern_cdf = function(x, R, veck) (1/(2*R))*as.numeric(x > -R)*(pmin(x ,R)+R)
    veck = 1
  } else {
    kk = degree/2-2
    area_row = vapply(0:(kk+2), FUN = function(i) 2*R^(2*i+1)/(2*i+1), FUN.VALUE = 1)
    zero_row = vapply(0:(kk+2), FUN = function(i) R^(2*i), FUN.VALUE = 1)
    deriv_row = c(0,vapply(0:(kk+1), FUN = function(i) 2*(i + 1)*R^(2*i+1), FUN.VALUE = 1))
    if (kk>0) {
      orth_rows = lapply(seq(0,max((2*kk-2),0),2), FUN = function(r) {
        vapply(0:(kk+2), FUN = function(i) 2*R^(2*i+3+r)/(2*i+3+r), FUN.VALUE = 1)
      })
      orth_rows = do.call(rbind, orth_rows) 
      mm = rbind(area_row, zero_row, deriv_row, orth_rows)
    } else mm = rbind(area_row, zero_row, deriv_row)
    
    mm_inv = solve(mm)
    veck = mm_inv %*% c(1, rep(0,kk+2))
    kern = function(x, R, veck) {
      ll = lapply(1:length(veck), FUN = function(c) veck[c]*x^(2*c-2))
      w = Reduce("+", ll)*(x > -R & x < R)
      return(w)
    }
    
    kern_cdf = function(x, R, veck) {
      u = pmin(x, R)
      ll = lapply(1:length(veck), FUN = function(c) veck[c]*(u^(2*c-1) + R^(2*c-1))/(2*c-1))
      w = Reduce("+", ll)*as.numeric(x > -R)
      return(w)
    }
  }
  
  return(list(veck = veck, R = R,  kern = kern, kern_cdf = kern_cdf))
}