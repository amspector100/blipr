#' A simulated point-source detection image from Liu et. al (2021).
#'
#' @format A list containing three attributes:
#' \describe{
#'   \item{image}{A 20 x 20 matrix representing an image.}
#'   \item{true_locs}{A 6 x 2 matrix representing the locations of sources.}
#'   \item{post_samples}{A 966 x 7 x 2 array of posterior samples from a StarNet model.}
#' }
#' @source \url{https://arxiv.org/pdf/2102.02409.pdf/}
"starnet_sim_data"
