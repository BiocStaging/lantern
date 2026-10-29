# SeqVarTools is imported only so that it installs with lantern:
# GMMAT::SMMAT() reads the DS dosages in our GDS files with
# SeqVarTools::imputedDosage(), but GMMAT only Suggests SeqVarTools.

#' lantern: Ancestry-Specific Rare Variant Analysis
#'
#' R package for ancestry-specific rare variant association analysis
#' with pure C backend for performance-critical operations.
#'
#' @section Step 1 -- split: \code{\link{ancestry_split}}
#' @section Step 2 -- write GDS: \code{\link{write_ancestry_gds}}
#' @section Step 3 -- associate: \code{\link{ancestry_smmat}}
#'
#' @keywords internal
#' @useDynLib lantern, .registration = TRUE
#' @importFrom utils head combn
#' @importFrom stats median pcauchy setNames
#' @importFrom SeqVarTools imputedDosage
"_PACKAGE"
