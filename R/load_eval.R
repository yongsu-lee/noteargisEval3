##########################################################
## Script for generating Evaluation Related Functions ####
##########################################################

skel = function(A){
  result_temp = A + t(A)
  A_skeleton = (result_temp != 0) * 1
  return(A_skeleton)
}


rev_edge = function(A_true, A_est){

  cpdag_A = pcalg::dag2cpdag(A_true)*1
  cpdag_est_A = pcalg::dag2cpdag(A_est)*1

  # find out undirected edges from true/estiamte graphs
  comp1 = cpdag_A == t(cpdag_A)
  comp2 = cpdag_est_A == t(cpdag_est_A)

  rev_edges_temp = ((A_est==1) & (t(A_true)==1))

  # Among the rev_edges_temp, we will exclude case that
  # ... corresponding edge is undirected in both
  # ... CPDAGs of estimated and true graph
  non_R = sum((rev_edges_temp & comp1) & comp2)
  return(sum(rev_edges_temp) - non_R)

}


#' Calculate Evaluation Metrics
#'
#' @param A_est_by_lam blah
#' @param A_true blah
#' @param n_edges blah
#'
#' @return blah
#' @export
#'
eval_by_lam = function(A_est_by_lam, A_true, n_edges){

  P = R = E = M = FP = SHD = JI = c()

  n_lams = length(A_est_by_lam)

  for (ell in 1:n_lams){ # ell = 2

    A_est = A_est_by_lam[[ell]]

    P[ell] = sum(A_est)
    if (ell > 1 & P[ell] == 0){P[ell] <- 1e-8}
    R[ell] = rev_edge(A_true, A_est)

    skel_A_true = skel(A_true)
    E[ell] = sum((skel_A_true == 1) & (A_est == 1)) - R[ell]

    if (ell > 1 & E[ell] == 0){E[ell] <- 1e-8}

    M[ell] = n_edges - E[ell] - R[ell]
    FP[ell] = P[ell] - R[ell] - E[ell]

    SHD[ell] = R[ell] + M[ell] + FP[ell]
    JI[ell] = E[ell] / (P[ell] + n_edges - E[ell])

  }

  func_result = list(P = P, R = R, E = E, M = M, FP = FP, SHD = SHD, JI = JI)
  return(func_result)

}

#' Select Model using the result from eval_by_lam
#'
#' @param one_simu_result blah
#' @param n_edges blah
#' @param sel_crit blah
#'
#' @return blah blah blah
#' @export
#'
mod_sel = function(one_simu_result, n_edges, sel_crit = "JI"){

  sel_idx = ifelse (sel_crit == "JI",  which.max(one_simu_result$JI),
                    which.min(one_simu_result$SHD))

  P= one_simu_result$P[sel_idx]
  E = one_simu_result$E[sel_idx]
  R = one_simu_result$R[sel_idx]
  # M = one_simu_result$M[sel_idx]
  FP = one_simu_result$FP[sel_idx]
  TPR = E/n_edges
  FDR = (R + FP)/P
  SHD = one_simu_result$SHD[sel_idx]
  JI = one_simu_result$JI[sel_idx]

  func_result = data.frame(P, E, R, FP, TPR, FDR, SHD, JI)

  return(func_result)

}
