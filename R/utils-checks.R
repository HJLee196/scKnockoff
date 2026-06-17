.check_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L ||
      is.na(x) || x <= 0 || x != as.integer(x)) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }
}

.check_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L ||
      is.na(x) || x < 0 || x != as.integer(x)) {
    stop("`", name, "` must be a nonnegative integer.", call. = FALSE)
  }
}

.check_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L ||
      is.na(x) || x <= 0) {
    stop("`", name, "` must be a positive number.", call. = FALSE)
  }
}

.check_nonnegative_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L ||
      is.na(x) || x < 0) {
    stop("`", name, "` must be a nonnegative number.", call. = FALSE)
  }
}

.check_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be either TRUE or FALSE.", call. = FALSE)
  }
}