npixels <- 20
post_samples <- scan("post_samples.txt")
post_samples <- aperm(array(post_samples, c(2, 7, 966)))
image <- scan("sim_image.txt")
image <- aperm(array(image, c(npixels, npixels, 2)))
true_locs <- scan("sim_locs.txt")
true_locs <- t(matrix(true_locs, 2, 6))
starnet_sim_data <- list(
  post_samples=post_samples,
  image=image,
  true_locs=true_locs
)
usethis::use_data(starnet_sim_data, overwrite=T)
