#' Forstmann et al.'s data
#'
#' A dataset containing the speed or accuracy manipulation for a Random Dot
#' Motion experiment.
#'
#' Details on the dataset can be found in the following paper:
#'
#' \strong{Striatum and pre-SMA facilitate decision-making under time pressure}
#'
#' Birte U. Forstmann, Gilles Dutilh, Scott Brown, Jane Neumann,
#' D. Yves von Cramon, K. Richard Ridderinkhof, Eric-Jan Wagenmakers.
#'
#' \emph{Proceedings of the National Academy of Sciences Nov 2008, 105 (45)
#' 17538-17542; DOI: 10.1073/pnas.0805903105}
#'
#' @format A data frame with 15818 rows and 5 variables:
#' \describe{
#'   \item{subject}{integer ID for each subject}
#'   \item{rt}{reaction time for each trial as a double}
#'   \item{condition}{Factor with 3 levels for Speed, Accuracy and
#'     Neutral}
#'   \item{stim}{Factor with 2 levels for Left and Right trials}
#'   \item{resp}{Factor with 2 levels for Left and Right responses}
#' }
#' @source \url{https://www.pnas.org/content/105/45/17538}
"forstmann"

#' A sampled object of a model of the Forstmann dataset
#'
#' A pmwgs object with a limited number of samples of the Forstmann dataset.
#'
#' The pmwgs object is missing one aspect, the pmwgs$data element. In order
#' to fully replicate the full object (ie to run more sampling stages) you will
#' need to add the data back in, via sampled_forstmann$data <- forstmann
#'
#' @section Samples Element:
#'
#' The samples element of a PMwG object contains the different types of samples
#' estimated by PMwG. These include the three main types of samples
#' \code{theta_mu}, \code{theta_sig} and \code{alpha} as well as a number of
#' other items which are detailed here.
#' \describe{
#'   \item{theta_mu}{samples used for estimating the model parameters (group
#'     level), an array of size (n_pars x n_samples)}
#'   \item{theta_sig}{samples used for estimating the parameter covariance
#'     matrix, an array of size (n_pars x n_pars x n_samples)}
#'   \item{alpha}{samples used for estimating the subject random effects, an
#'     array of size (n_pars x n_subjects x n_samples)}
#'   \item{stage}{A vector containing what PMwG stage each sample was drawn in}
#'   \item{subj_ll}{The winning particles log-likelihood for each subject and
#'     sample}
#'   \item{a_half}{Mixing weights used during the Gibbs step when creating a
#'     new sample for the covariance matrix}
#'   \item{last_theta_sig_inv}{The inverse of the last samples covariance
#'     matrix}
#'   \item{idx}{The index of the last sample drawn}
#' }
#'
#' @format A pmwgs object minus the data. A pmwgs opbject is a list with a
#' specific structure and elements, as outlined below.
#' \describe{
#'   \item{par_names}{A character vector containing the model parameter names}
#'   \item{n_pars}{The number of parameters in the model}
#'   \item{n_subjects}{The number of unique subject ID's in the data}
#'   \item{subjects}{A vector containing the unique subject ID's}
#'   \item{prior}{A list that holds the prior for \code{theta_mu} (the model
#'     parameters). Contains the mean (\code{theta_mu_mean}), covariance matrix
#'     (\code{theta_mu_var}) and inverse covariance matrix
#'     (\code{theta_mu_invar})}
#'   \item{ll_func}{The log likielihood function used by pmwg for model
#'     estimation}
#'   \item{samples}{A list with defined structure containing the samples, see
#'     the Samples Element section for more detail}
#' }
#' @source \url{https://www.pnas.org/content/105/45/17538}
"sampled_forstmann"
