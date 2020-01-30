

#' Compute a single timestep of a linear time invariant system with white noise disturbances and measurement error.
#'
#' @param x x_t, a n length vector of observations
#' @param A An n by n matrix of system coefficients 
#' @param B An n by matrix of control-to-node mapping coefficients
#' @param u An m length vector of control inputs
#' @param theta A positive definite covariance matrix for white noise disturbances
#' @param gamma A positive definite covariance matrix for white noise measurement error
#'
#' @return A list with x: a vector of the states at t+1, and y: a vector of the observations at t+1
#' @export
#'
#' @examples
control_step <- function(x, A, B, u, theta = NA, gamma = NA){

  x_t1 = A %*% x + B %*% u

  if(is.matrix(theta)){

    x_t1 = x_t1 + mvrnorm(mu = rep(0, length(x)), Sigma = theta)
  }

  if(is.matrix(gamma)){

    y_t1 = mvrnorm(mu = x_t1, Sigma = gamma)
    return(list(x = x_t1, y = y_t1))
  }else{
    return(list(x = x_t1, y = x_t1))

  }

}

#' Simulate the trajectory of a linear time invariant system with white noise disturbances and measurement error, subject to control inputs.
#'
#' @param t_max The number of time points to simulate.
#' @param x_0 The starting values for the state vector.
#' @param B An n by matrix of control-to-node mapping coefficients
#' @param u An m length vector of control inputs
#' @param theta A positive definite covariance matrix for white noise disturbances
#' @param gamma A positive definite covariance matrix for white noise measurement error
#' @param G_seq A single Kalman gain matrix, or a sequence of Kalman gain matrices.
#' @param J_func A quadratic regulator cost function.
#' @param S The S matrix.
#' @param Q_seq A single Q matrix or a sequence of Q matrices
#' @param R_seq A single R matrix of a sequence of R matrices
#' @param delta For saturated inputs, a length 2 vector of the threshold at which the input is saturated, and the value of a saturated input. NA for no saturated inputs
#' @param d_nosign Logical, if inputs are saturated, should the saturation only occur for positive inputs.
#' @param d_toggle Logical, if inputs are saturated, should the input be a step function. (0 unless it exceeds the saturation point.)
#'
#' @return A list with state values, observation values, control input values and cost function values for 1:t_max timepoints.
#' @export
#'
#' @examples
control_traj <- function(t_max, x_0, A, B, theta, gamma, G_seq, J_func, S, Q_seq, R_seq, delta = NA, d_nosign = F, d_toggle = F){


  if(any(!(c(length(G_seq),length(Q_seq),length(R_seq))%in% c(1, t_max)))){
    simpleError("All matrix sequences must be 1 or t_max in length. Sequences can be differently 1 or t_max in length")
  }

  j_vec = vector(length =t_max-1)

  x_mat = matrix(0, t_max, length(x_0))
  x_mat[1,] = x_0
  y_mat = matrix(0, t_max, length(x_0))

  if(!is.na(gamma)){

    y_mat[1,] = mvrnorm(mu = x_0, Sigma = gamma)
  }else{
    y_mat[1,] = x_0
  }


  S_seq = S_seq_calc(A, S, B, Q_seq, R_seq, t_max)
  u_mat = matrix(0,t_max-1 ,dim(B)[[2]])
  for(i in 1:(t_max-1)){
    if(is.matrix(G_seq)){G = G_seq}else{G = G_seq[[i]]}
    u = G %*% y_mat[i,]
    print(u)
    if(!is.na(delta)){
      u = sign(u)*(abs(u) >= delta[1])*delta[2] + u*(abs(u) < delta[1])

      if(d_toggle){
        u =sign(u)*(abs(u) >= delta[1])*delta[2]
      }
      if(d_nosign){
        u = (u > delta[1])*delta[2]+ u*(u <= delta[1])*(u > 0)
      }


    }

    temp = control_step(x_mat[i,], A, B, u, theta, gamma)
    x_mat[i+1,] = t(temp$x)
    y_mat[i+1,] = t(temp$y)
    u_mat[i,] = t(u)
    j_vec[i] = perf_index(x_mat[i+1,],S_seq[[i]])
  }

  return(list(x_mat,y_mat,u_mat,j_vec))

}

#' Recursively solve the Riccati equations
#'
#' @param A 
#' @param S 
#' @param B 
#' @param Q_seq 
#' @param R_seq 
#' @param tmax 
#'
#' @return
#' @export
#'
#' @examples
S_seq_calc <- function(A, S, B, Q_seq, R_seq, tmax){

  S_seq = list()
  S_seq[[tmax]] = S

  for(i in (tmax-1):1){

    if(is.matrix(Q_seq)){Q = Q_seq}else{Q = Q_seq[[i]]}
    if(is.matrix(R_seq)){R = R_seq}else{R = R_seq[[i]]}
    S_seq[[i]] = t(A)%*%(S_seq[[i+1]] - S_seq[[i+1]] %*% B %*% solve(t(B) %*% S_seq[[i+1]] %*% B + R) %*% t(B) %*% S_seq[[i+1]]) %*% A + Q
  }

  return(S_seq)

}

perf_index <- function(x, S_i){

  j_i = .5*t(x) %*% S_i %*% x

  return(j_i)

}

G_seq_calc <- function(S_seq, R_seq, A, B, t_max){
  
  G_seq = list()
  for(i in (t_max-1):1){
    if(is.matrix(R_seq)){R = R_seq}else{R = R_seq[[i]]}
    G_seq[[i]] = -solve(t(B) %*% S_seq[[i+1]]%*% B + R) %*% t(B) %*%S_seq[[i+1]] %*% A
  }
  return(G_seq)

}

optimal_control_closed_loop <- function(t_max, A, B, S, Q_seq, R_seq){

  S_seq = S_seq_calc(A, S, B, Q_seq, R_seq, t_max)
  G_seq = G_seq_calc(S_seq, R_seq, A, B, t_max)

  return(G_seq)
}



