#' @title Reynolds Adolescent Depression Scale (RADS-2) responses of Peruvian adolescents
#'
#' @description Responses of 917 adolescents from Lima (Peru) to the 25 items retained in the four-factor version of the
#' Reynolds Adolescent Depression Scale, second edition (RADS-2), each answered on a four-point Likert scale. The data come
#' from Ramos-Vera et al. (2023), who also grouped the items into four clusters (stored with the data).
#'
#' @format A data frame with 917 rows (adolescents) and 25 integer columns (items), with values 1 to 4. Higher values
#' indicate more frequent depressive symptoms for every item. The columns keep the item numbering of the RADS-2
#' (\code{D1}, \code{D3}, ..., \code{D30}). The attribute \code{"clusters"} is a named character vector giving the
#' cluster of each item:
#' \describe{
#'   \item{Dysphoria}{D3, D6, D7, D8, D16, D21, D26}
#'   \item{Anhedonia/Negative Affect}{D1, D5, D12}
#'   \item{Negative Self-Assessment}{D4, D9, D13, D14, D15, D19, D20, D30}
#'   \item{Somatic Complaint}{D11, D17, D18, D22, D24, D27, D28}
#' }
#' Subsetting the data frame (e.g. \code{rads2[, 1:10]}) drops the attribute; retrieve it first with
#' \code{attr(rads2, "clusters")}.
#'
#' @details The data are derived from "S1 Data" of Ramos-Vera et al. (2023), with the following changes: the 25 items of
#' the four-factor version were kept (items 2, 10, 23, 25 and 29, and the age and sex of the respondents, were removed),
#' the items were renamed from \code{A<k>} to \code{D<k>}, and the positively worded items \code{D1}, \code{D5} and
#' \code{D12} were reverse-coded (\code{5 - x}), so that higher values indicate more symptoms for every item. The rows
#' are in the original order. The item wording of the RADS-2 is not included.
#'
#' @source Ramos-Vera, C., Quispe Callo, G., Basauri Delgado, M., Vallejos Saldarriaga, J., and Saintila, J. (2023).
#' "S1 Data. Reynolds Adolescent Depression Scale", supporting information of the article below,
#' \doi{10.1371/journal.pone.0286081.s001}. Licensed under the Creative Commons Attribution 4.0 International license
#' (CC BY 4.0, \url{https://creativecommons.org/licenses/by/4.0/}); the changes are listed under Details. The rest of
#' the package is licensed under the MIT license.
#'
#' @references Ramos-Vera, C., Quispe Callo, G., Basauri Delgado, M., Vallejos Saldarriaga, J., and Saintila, J. (2023).
#' Factorial and network structure of the Reynolds Adolescent Depression Scale (RADS-2) in Peruvian adolescents.
#' \emph{PLOS ONE}, 18(5), e0286081. \doi{10.1371/journal.pone.0286081}
#'
#' Reynolds, W. M. (2002). \emph{Reynolds Adolescent Depression Scale, 2nd Edition: Professional Manual}. Lutz, FL:
#' Psychological Assessment Resources.
#'
#' @examples
#' data(rads2)
#' dim(rads2)
#' table(attr(rads2, "clusters"))
#'
#' # point estimates and sandwich standard errors for the seven dysphoria items
#' dysphoria <- names(which(attr(rads2, "clusters") == "Dysphoria"))
#' fit <- dmrfit(rads2[, dysphoria], with_prior = TRUE)
#' summary(fit)
"rads2"
