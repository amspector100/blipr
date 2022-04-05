inverse_size <- function(cand_group) {
	return(1 / length(cand_group$group))
}

log_inverse_size <- function(cand_group) {
	return(1 + log(inverse_size(cand_group), 2))
}

inverse_radius <- function(cand_group) {
  return(1 / cand_group$radius)
}
