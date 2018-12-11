init_sums <- function(x_train)

update_sums <- function(sums, new.x, t, w) {
  t2 <- min(w + 4, t)
  previous.w <- sums$w[, t2 - 2]
  previous.v <- sums$v[, t2 - 2]
  new.w <- previous.w + new.x
  new.deviation <- new.x - 1 / (t - 1) * previous.w
  new.v <- previous.v + (t - 1) / t * new.deviation^2
  new.s <- new.v / t
  if (t <= w + 3) {
    sums$w[, t2 - 1] <- new.w
    sums$v[, t2 - 1] <- new.v
    sums$s[, t2 - 1] <- new.s
  } else {
    sums$w <- cbind(sums$w[, 2:(t2 - 2), drop = FALSE], new.w)
    sums$v <- cbind(sums$v[, 2:(t2 - 2), drop = FALSE], new.v)
    sums$s <- cbind(sums$s[, 2:(t2 - 2), drop = FALSE], new.s)
  }
  return(sums)
}

bartlett_corr <- function(t, w) {
  k <- max(2, t - w - 1):(t - 2)
  1 / 2 * (-t * log(t) + k * log(k) + (t - k) * log(t - k) +
           t * digamma((t - 1) / 2) - k * digamma((k - 1) / 2) -
           (t - k) * digamma((t - k - 1) / 2))
}

mixture_log_liks <- function(sums, t, w, p0, n.liks = min(t - 3, w)) {
  t2 <- min(w + 3, t)
  w.full <- sums$w[, t2 - 1]
  v.full <- sums$v[, t2 - 1]
  s.full <- sums$s[, t2 - 1]
  N <- length(w.full)
  ks <- max(2, t - (w + 1)):(t - 2)
  log.likelihoods <- matrix(0, nrow = N, ncol = n.liks)
  for (i in 1:n.liks) {
    j <- max(t2 - 3 - n.liks, 0) + i
    k <- ks[j]

    s.pre <- sums$s[, j]

    w.pre <- sums$w[, j]
    avg.pre <- w.pre / k
    avg.post <- (w.full - w.pre) / (t - k)
    v.pre <- sums$v[, j]
    v.post <- v.full - v.pre - (k * (t - k)) / t * (avg.pre - avg.post)^2
    s.post <- 1 / (t - k) * v.post
    if (any(s.post <= 0)) {
      warning(paste('Numerical error in estimating post-change variance', t, t - k))
      s.post[which(s.post <= 0)] <- NA
    }
    log.likelihoods[, i] <- t * log(s.full) - k * log(s.pre) - (t - k) * log(s.post)
  }

  corrected.log.liks <- matrix(0, nrow = N, ncol = n.liks)
  C <- CalcBartlettVar(t, w)
  for (j in 1:N) {
    corrected.log.liks[j, ] <- log.likelihoods[j, ] / C
  }
  global.log.liks <- apply(corrected.log.liks, 2,
                           function(x) {sum(log(1 - p0 + p0 * exp(1/2 * x)))})
  return(global.log.liks)
}
