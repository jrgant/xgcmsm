
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{age} and \code{sqrt.age}. The \code{sqrt.age} vertex attribute is also
#' updated on the three networks.
#'
#' @keywords module msm
#' @export
#'
aging_msm <- function(dat, at) {

  age.wk <- dat$attr$age.wk
  age <- dat$attr$age
  age.grp <- dat$attr$age.grp
  active <- dat$attr$active

  age.wk[active == 1] <- age.wk[active == 1] + 1
  age[active == 1] <- age.wk[active == 1] / 52

  age.grp[active == 1] <- cut(
    age[active == 1],
    dat$param$netstats$demog$age.breaks,
    right = FALSE,
    labels = FALSE
  )

  dat$attr$age.wk <- age.wk
  dat$attr$age.grp <- age.grp
  dat$attr$age <- age

  return(dat)
}
