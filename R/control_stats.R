#' Controllability Gramian
#'
#' Compute the (infinite time) controllability Gramian for the discrete linear time invariant system described by \eqn{x(t+1) = Ax(t) + Bu(t)}. 
#' The infinite time controllability Gramian is the solution to the discrete Lyapunov equation \eqn{AWA^\prime-W = -BB^\prime}, while the finite time Gramian for time \eqn{T} is 
#' \deqn{W_t = \sum_{t = 0}^T A^tBB^\prime(A^\prime)^t} 
#'
#' @param A \eqn{n x n} matrix.
#' @param B \eqn{n x m} matrix.
#' @param t Either NA for infinite time Gramian, or a positive non-zero integer. Defaults to NA.
#' 
#' @return The infinite time or finite time controllability Gramian
#' @export
#'
#' @examples
#' 
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' B = diag(3)
#' 
#' #Infinite time Gramian
#' W_inf = control_gramian(A, B)
#' 
#' #4 time Gramian
#' W_4 = control_gramian(A,B,4)
control_gramian <-function(A,B, t = NA){
  
  W = B%*%t(B)
  if(is.na(t)){
  gramian = dlyap(A, W)
  }else{
  gramian = matrix(0, dim(A)[[1]], dim(A)[[2]])
  
  for(i in 1:t){
    
    gramian = gramian = A%^%t %*% (-W) %*% (t(A))%^%t
    
  }
  }
  return(gramian)
}



#' Trace of the Inverse Gramian
#' 
#' A commonly used measure of the overall controllability of a system defined by \eqn{x_(t+1) = Ax_(t) + Bu_(t)}.
#'
#' @param A An n by n matrix.
#' @param B An n by m matrix.
#'
#' @return Trace of the inverse infinite time Gramian.
#' @export
#'
#' @examples
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' B = diag(3)
#' 
#' inv_average_control(A, B)
inv_average_control <- function(A,B){
  
  gramian <- control_gramian(A, B)
  return(sum(diag(pinv(gramian))))
  
}

#' Trace of the  Gramian
#' 
#' A commonly used measure of the overall controllability of a system defined by \eqn{x_(t+1) = Ax_(t) + Bu_(t)}.
#'
#' @param A An n by n matrix.
#' @param B An n by m matrix.
#'
#' @return Trace of the infinite time Gramian.
#' @export
#'
#' @examples
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' B = diag(3)
#' 
#' average_control(A, B)
average_control <- function(A, B){
  gramian <- control_gramian(A, B)
  return(sum(diag(gramian)))
  
}

#' Modal Control
#' 
#' Calculates the modal control of a system defined by \eqn{x_(t+1) = Ax_(t) + Bu_(t)}.
#'
#' @param A An n by n matrix.
#' @param B An n by m matrix.
#'
#' @return An m by n matrix representing the control of the nth mode by the mth control input.
#' @export
#'
#' @examples
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' B = diag(3)
#' 
#' modal_control(A, B)
modal_control <- function(A, B){
  
  
  eig_vector = eigen(A/(1+svd(A)$d[1]))$vector
  
  if(is.vector(B)){
    B = matrix(B, length(B), 1)
    
  }
  
  output = matrix(0, dim(A)[[1]], dim(B)[[2]])
  
  for(i in 1:dim(eig_vector)[[1]]){
    for(j in 1:dim(B)[[2]]){
      
      output[i,j] = Mod(t(eig_vector[,i])%*%B[,j])/(sqrt(sum(Mod(eig_vector[,i])^2))*sqrt(sum(B[,j]^2)))
      
      
    }
  }
  return(output)
}

#' Modal Control Centrality
#' 
#' Calculates the modal control centrality of a system defined by \eqn{x_(t+1) = Ax_(t)}.
#'
#' @param A An n by n matrix.
#'
#' @return A length n vector of modal control centrality measures, representing the overall modal control of each node in the system.
#' @export
#'
#' @examples
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' 
#' modal_control_centrality(A, B)
modal_control_centrality <- function(A){
  
  eig = eigen(A/(1+svd(A)$d[1]))
  eig_vector = eig$vector
  eig_vals = eig$values
  
  
  output = matrix(0, dim(A)[[1]], dim(A)[[2]])
  
  for(i in 1:dim(eig_vector)[[1]]){
    for(j in 1:dim(eig_vector)[[1]]){
      
      output[j,i] = Mod(t(eig_vector[i,j]))^2*(1-Mod(eig_vals[j])^2)
      
      
    }
  }
  return(colSums(output))
}

#' Average Control Centrality
#' 
#' Calculates the modal control centrality of a system defined by \eqn{x_(t+1) = Ax_(t) + Bu_(t)}.
#'
#' @param A An n by n matrix.
#'
#' @return A length n vector of average control centrality measures, representing the overall average control of each node in the system.
#' @export
#'
#' @examples
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' 
#' ave_control_centrality(A)
ave_control_centrality <- function(A){
  
  centrality = diag(control_gramian(t(A)/(1+Mod(eigen(t(A))$values[1])), diag(dim(A)[[1]])))
  
  return(centrality)
}