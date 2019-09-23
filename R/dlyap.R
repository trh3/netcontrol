library(Matrix)
library(pracma)

#
#A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
# 
# C = matrix(c(-2,-8,11,2,-6,13,-3,-5,-2), 3,3)
# 
# X = dlyap(t(A), C)
# 
# 
# 
# W = diag(1, 3)
#
#
#
#

#' Discrete Lyapunov Equation Solver, base R implementation
#'
#' @param A 
#' @param W 
#'
#' @return
.dlyap_native <- function(A, W){
  
  
  Sch_A <-Schur(t(A), vectors = T)
  
  # From Datta 2004 pg. 275, "Numerical Solutions and Conditioning of Lyapunov and Sylvester Equations"
  X_prime = matrix(0, dim(A)[[1]], dim(A)[[2]])
  Sch_T = Sch_A$T
  Q = Sch_A$Q
  C = t(Q)%*%W%*%Q
  
  i = dim(A)[[1]]
  while(i != 0){
  
    if(i == dim(A)[[1]]){
      if(Sch_T[i,i-1] == 0){
        X_prime[,i] = solve(Sch_T[i,i] * Sch_T - diag(1, dim(A)[[1]])) %*% C[,i]
        i = i-1
      }else{
        temp = solve((Sch_T[(i-1):i,(i-1):i] %x% Sch_T) -diag(1, dim(A)[[1]]*2))
        c = as.numeric(C[,((i-1):i)])
        X_prime_temp = as.matrix(temp %*% c, dim(A)[[1]], 2)
        X_prime[,(i-1):i] = X_prime_temp
        i = i-2
      }
    }else{
      if(i==1){
        temp = solve(Sch_T[i,i] * Sch_T - diag(1, dim(A)[[1]]))
        sub = 0
        for(j in (i+1): dim(A)[[1]]){
          sub = sub + Sch_T[i,j]*X_prime[,j]
        }
        sub = Sch_T%*% sub
        
        
        sub = C[,i] - sub
        X_prime[,i] = temp %*% sub
        i = i -1
      }else{
        if(Sch_T[i,i-1] == 0){
          temp = solve(Sch_T[i,i] * Sch_T - diag(1, dim(A)[[1]]))
          sub = 0
          for(j in (i+1): dim(A)[[1]]){
            sub = sub + Sch_T[i,j]*X_prime[,j]
          }
          sub = Sch_T%*% sub
          
          
          sub = C[,i] - sub
          X_prime[,i] = temp %*% sub
          i = i-1
        }else{
          temp = solve((Sch_T[(i-1):i,(i-1):i] %x% Sch_T) -diag(1, dim(A)[[1]]*2))
          sub_1 = 0
          sub_2 = 0
          for(j in (i+1): dim(A)[[1]]){
            sub_1 = sub_1 +Sch_T[i,j]*X_prime[,j]
            sub_2 = sub_2 +Sch_T[i-1,j]*X_prime[,j]
          }
          sub_1 = Sch_T%*% sub_1
          sub_2 = Sch_T%*% sub_2
          
          c_1 = as.numeric(C[,(i)]) - sub_1
          c_2 = as.numeric(C[,(i-1)]) - sub_2
          c = rbind(c_2,c_1)
          
          X_prime_temp = as.matrix(temp %*% c, dim(A)[[1]], 2)
          X_prime[,(i-1):i] = X_prime_temp
          i = i-2
        }
        
        
      }
      
    }
    
    
  }
  
  X = Q%*%X_prime%*%t(Q)
  
  return(X)
  
}

#' Discrete Lyapunov Equation Solver
#' 
#' Computes the solution of \eqn{AXA^T - X + W = 0} using the Barraud 1977 approach, adapted from Datta 2004. This implementation is equivalent to the Matlab implementation of dylap.
#'
#' @param A \eqn{n x n} numerical or complex matrix.
#' @param W \eqn{n x n} numerical or complex matrix.
#'
#' @return The solution to the above Lyapunov equation.
#' @export
#'
#' @references
#' \insertRef{barraud_numerical_1977}{netcontrol}
#' 
#' \insertRef{datta_numerical_2004}{netcontrol}
#' 
#' @examples 
#' A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
#' C = matrix(c(-2,-8,11,2,-6,13,-3,-5,-2), 3,3)
#' X = dlyap(t(A), C)
#' 
#' print(sum(abs(A %*% X %*% t(A) - X + C)))
dlyap <- function(A, W){
  
  W = -W
  Sch_A <-Schur(t(A), vectors = T)
  
  block_diag = !is.complex(eigen(t(A))$values)
  # From Datta 2004 pg. 275, "Numerical Solutions and Conditioning of Lyapunov and Sylvester Equations"
  X_prime = matrix(0, dim(A)[[1]], dim(A)[[2]])
  Sch_T = Sch_A$T
  Q = Sch_A$Q
  C = t(Q)%*%W%*%Q
  X = dlyap_internal(Sch_T, Q, C, dim(A)[[1]], block_diag)
  return(X)
  
}
