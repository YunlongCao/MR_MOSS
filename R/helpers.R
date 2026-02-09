mrmoss_null_if_empty <- function(x) {
  if (length(x) == 0 || all(is.na(x)) || all(trimws(as.character(x)) == "")) {
    return(NULL)
  }
  x
}

mrmoss_to_numeric <- function(x, default = NA_real_) {
  x <- mrmoss_null_if_empty(x)
  if (is.null(x)) {
    return(default)
  }
  suppressWarnings(as.numeric(x)[1])
}

mrmoss_to_character <- function(x, default = NULL) {
  x <- mrmoss_null_if_empty(x)
  if (is.null(x)) {
    return(default)
  }
  as.character(x)[1]
}

mrmoss_expand_env_vars <- function(path) {
  path <- as.character(path)[1]
  if (!grepl("\\$\\{", path)) {
    return(path)
  }

  pattern <- "\\$\\{([A-Za-z0-9_]+)\\}"
  vars <- unique(gsub(pattern, "\\1", regmatches(path, gregexpr(pattern, path))[[1]]))

  for (v in vars) {
    token <- paste0("${", v, "}")
    value <- Sys.getenv(v, unset = "")
    if (identical(value, "")) {
      stop(sprintf("Environment variable %s is not set but required by path: %s", v, path), call. = FALSE)
    }
    path <- gsub(token, value, path, fixed = TRUE)
  }

  path
}

mrmoss_resolve_path <- function(path, root_dir = NULL) {
  path <- mrmoss_expand_env_vars(path)

  if (!is.null(root_dir) && !grepl("^/", path)) {
    path <- file.path(root_dir, path)
  }

  normalizePath(path, winslash = "/", mustWork = FALSE)
}

mrmoss_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

mrmoss_get_plink_binary <- function(plink_bin = NULL) {
  if (!is.null(plink_bin) && nzchar(plink_bin)) {
    return(plink_bin)
  }

  if (requireNamespace("plinkbinr", quietly = TRUE)) {
    ns <- asNamespace("plinkbinr")
    if (exists("get_plink_exe", envir = ns, mode = "function")) {
      return(get("get_plink_exe", envir = ns)())
    }
    if (exists("get_plink_binary", envir = ns, mode = "function")) {
      return(get("get_plink_binary", envir = ns)())
    }
  }

  stop(
    "No PLINK binary provided. Set `plink_bin` or install `plinkbinr`.",
    call. = FALSE
  )
}

mrmoss_make_spd <- function(R, eps = 1e-6) {
  R <- (R + t(R)) / 2
  diag(R) <- 1

  eig <- eigen(R, symmetric = TRUE)
  if (min(eig$values) >= eps) {
    return(R)
  }

  eig$values[eig$values < eps] <- eps
  R_spd <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  D_inv <- diag(1 / sqrt(diag(R_spd)), nrow = nrow(R_spd))
  R_spd <- D_inv %*% R_spd %*% D_inv
  R_spd <- (R_spd + t(R_spd)) / 2
  diag(R_spd) <- 1
  R_spd
}
