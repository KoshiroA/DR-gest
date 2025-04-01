

if (requireNamespace("doParallel", quietly = TRUE)) {
  suppressWarnings(library(doParallel))
} else {
  stop("Required package 'doParallel' is not installed.")
}

if (requireNamespace("doRNG", quietly = TRUE)) {
  suppressWarnings(library(doRNG))
} else {
  stop("Required package 'doRNG' is not installed.")
}

if (requireNamespace("brm", quietly = TRUE)) {
  suppressWarnings(library(brm))
} else {
  stop("Required package 'brm' is not installed.")
}

if (requireNamespace("geeM", quietly = TRUE)) {
  suppressWarnings(library(geeM))
} else {
  stop("Required package 'geeM' is not installed.")
}

if (requireNamespace("tidyverse", quietly = TRUE)) {
  suppressWarnings(library(tidyverse))
} else {
  stop("Required package 'tidyverse' is not installed.")
}
