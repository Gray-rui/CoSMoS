
mod_simulateTS <- function (aTS, from = NULL, to = NULL) {
  dist <- attr(aTS, "dist")
  acsID <- attr(aTS, "acsID")
  season <- attr(aTS, "season")
  date <- data.table(attr(aTS, "date"))
  x <- aTS$data
  r <- reportTS(aTS, method = "stat")
  f <- aTS$dfits
  a <- aTS$afits
  ACS <- list()
  if (dist == "beta") {
    distbounds = c(0, 1)
  } else {
    distbounds = c(-Inf, Inf)
  }
  for (i in seq_along(x[[1]])) {
    p <- actpnts(margdist = attr(f[[i]], "dist"), margarg = f[[i]], 
                 p0 = r[i, "p0"], distbounds = distbounds)
    fp <- fitactf(p)
    lag <- 0:(length(attr(a[[i]], "eACS")) - 1)
    id <- attr(a[[i]], "ID")
    as <- do.call(acs, c(list(id = id, t = lag), a[[i]]))
    ACS[[i]] <- actf(as, fp$actfcoef[1], fp$actfcoef[2])
  }
  names(ACS) <- names(a)
  p0 <- uval <- gauss <- value <- . <- season_id <- rn <- NULL
  if (is.null(from)) {
    from <- date[1, date]
  }
  if (is.null(to)) {
    to <- date[.N, date]
  }
  by <- difftime(date[2, date], date[1, date])
  gausian <- mod_seasonalAR(x = seq(from = from, to = to, by = by), 
                            ACS = ACS,
                            season = season)
  setkey(x = gausian, season)
  para <- as.data.table(x = t(sapply(f, function(x) {
    as.matrix(x = do.call(what = cbind, x))
  })))
  names(para) <- names(f[[1]])
  para[, `:=`(season, as.numeric(gsub("data_nz_", "", rownames(para))))]
  para[, `:=`(p0, r[, "p0"])]
  setkey(para, season)
  aux <- merge(gausian, para, all.x = T)
  aux <- aux[order(date)]
  aux[, `:=`(uval, (pnorm(q = gauss) - p0)/(1 - p0))]
  aux[uval < 0, `:=`(uval, 0)]
  d <- getDistArg(dist)
  for (i in para[, season]) {
    trans.para <- para[season == i, !c("p0", "season")]
    aux[season == i, `:=`(value, do.call(what = paste0("q", 
                                                       dist), args = c(list(p = uval), as.list(trans.para))))]
  }
  out <- aux[, .(date, value)]
  return(out)
}



mod_seasonalAR <- function (x, ACS, season = "month") {
  time <- data.table(time = x)
  y <- s <- n <- . <- id <- value <- NULL
  time[, `:=`(y, year(time))]
  time[, `:=`(s, do.call(season, .(time)))]
  time[, `:=`(n, .N), by = .(y, s)]
  d <- as.data.frame(unique(time[, -1]))
  alpha <- lapply(ACS, YW)
  out <- data.table(value = AR1(max(sapply(alpha, length)), 
                                ACS[[which(gsub(paste(season, ""), "", names(ACS)) == 
                                             d[1, "s"])]][2]))
  out[, `:=`(id, 0)]
  esd <- c()
  for (i in seq_along(alpha)) {
    esd[i] <- sqrt(1 - sum(alpha[[i]] * ACS[[i]][2:length(ACS[[i]])]))
  }
  for (j in 1:dim(d)[1]) {
    s <- d[j, "s"]
    ss <- which(gsub(paste(season, ""), "", names(alpha)) == 
                  s)
    p <- length(alpha[[ss]])
    val <- out[(dim(out)[1] + 1 - p):dim(out)[1], value]
    aux <- length(val)
    n <- unlist(d[j, "n"])
    if (season == "year"){
      gn <- rnorm(n + p, mean = 0, sd = esd[j])
    } else {
      gn <- rnorm(n + p, mean = 0, sd = esd[s])
    }
    a.rev <- rev(alpha[[ss]])
    for (i in (p + 1):(n + p)) {
      val[i] <- sum(val[(i - p):(i - 1)] * a.rev) + gn[i]
    }
    out <- rbind(out, data.table(value = val[-1:-aux], id = s))
  }
  out <- out[id != 0, ]
  dt <- data.table(date = x, gauss = out[, value], season = out[, 
                                                                id])
  return(dt)
}
