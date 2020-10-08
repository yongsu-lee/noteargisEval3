##########################################################
## Script for generating cd algo related functions #######
##########################################################


#' Convert data frame as sparsebnData object
#'
#' @param data blah
#'
#' @return blah
#' @export
#' @importFrom sparsebnUtils sparsebnData
#'
conv_to_cd_data = function(data){
  func_result = sparsebnUtils::sparsebnData(data, type = "discrete")
  return(func_result)
}



#' Calculate an Adjacency Matrix using sparsebn object
#'
#' @param result_cd blah
#'
#' @return blah
#' @export
#' @importFrom sparsebnUtils get.adjacency.matrix get.lambdas
calc_A_est_by_lam_cd = function(result_cd){

  n_lams = length(sparsebnUtils::get.lambdas(result_cd))
  func_result = list(0)

  for (ell in 1:n_lams){ # ell = 1

    func_result[[ell]] = as.matrix(
      sparsebnUtils::get.adjacency.matrix(result_cd[[ell]]) )

  }

  return(func_result)

}
