
# von Mises ---------------------------------------------------------------

rnpvm = function(n=100, mu=0, pr=1, beta=1, mix, family=FALSE) {
  if (missing(mix)) mix = disc(mu, pr)
  if (n == 0) return(numeric(0))
  l = length(mix$pt)
  if(length(beta) != l) beta = rep(beta, length.out = l)
  pr = pr / sum(pr)
  suppressWarnings(index <- sample.int(l, n, prob = mix$pr, replace = TRUE))
  value = rep(0, n)
  for (i in 1:n) value[i] =
    circular::rvonmises(1, circular::circular(mix$pt[index[i]]), beta[index[i]])
  x = list(v = value, w = 1)
  class(x) = "npvm"
  if (family) x  
  else value
}

set.seed(17891642)
x1 = rnpvm(100, mu = c(pi, pi/2), pr = c(0.3, 0.7), beta = 20)

as.npvm = function(x, w=1) {
  x = list(v = c(x), w = w)
  class(x) = "npvm"
  x
}

npvm = function(v, w) {
  x = list(v = v, w = w)
  class(x) = "npvm"
  x
}

length.npvm = function(x) length(x$v)

weight.npvm = function(x, beta) x$w

gridpoints.npvm = function(x, beta, grid=100) {
  seq(0, 2 * base::pi, length = grid)
}

initial.npvm = function(x, beta, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = 20
  if(is.null(mix) || is.null(mix$pt)) {
    breaks = seq(0, 2 * base::pi, length = 20)
    r = nspmix::whist(x$v, x$w, breaks = breaks, freq = FALSE, plot = FALSE)
    i = r$density != 0
    mix = disc(r$mids[i], r$density[i])
  }
  list(beta = beta, mix = mix)
}

valid.npvm = function(x, beta, theta) beta > 0

suppspace.npvm = function(x, beta) c(0, 2 * base::pi)

logd.npvm = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","dk","dt")

  if(which[1] == 1){
    dl$ld = beta * cos(x$v - rep(pt, rep(length(x$v), length(pt)))) -
      log(2 * base::pi * besselI(beta, 0, expon.scaled = TRUE)) - beta
    dl$ld[dl$ld < -300] = -300  # if density is zero, change to 1e-100, log density is -230
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  if(which[2] == 1)
    dl$db = NULL
  if(which[3] == 1){
    dl$dt = beta * sin(x$v - rep(pt, rep(length(x$v), length(pt))))
    dim(dl$dt) = c(length(x$v), length(pt))
  }
  dl
}

dvm1 = function(x, mu, beta, log=FALSE) {
  lm = length(mu)
  lx = length(x)
  if(length(beta) != lm) beta = rep(beta, length.out = lm)
  x.v = rep(x, lm)
  mu.v = rep(mu, rep(lx, lm))
  beta.v = rep(beta, rep(lx, lm))
  b0 = besselI(unique(beta.v), 0, expon.scaled = TRUE)
  b0.v = rep(rep(b0, length.out = lm), rep(lx, lm))
  kcos = beta.v * cos(x.v - mu.v)
  constant = 2 * pi * b0.v
  if(log) d = kcos - log(constant) - beta.v
  else d = exp(kcos) / constant / exp(beta.v)
  dim(d) = c(lx, lm)
  d
}
dvm1(1, 2, 20)
dvm1(1:2, 2, 20)
dvm1(1, 2:3, 20)
dvm1(1:10, 1:2, 20)

dvm = function(x, mu, beta, log=FALSE) {
  lm = length(mu)
  lx = length(x)
  lb = length(beta)
  maxlen = max(lm, lx, lb)  # repeat to the largest length
  if(lx != maxlen) x = rep(x, length.out = maxlen)
  if(lm != maxlen) mu = rep(mu, length.out = maxlen)
  b0 = besselI(unique(beta), 0, expon.scaled = TRUE)
  if(lb != maxlen) b0 = rep(b0, length.out = maxlen)
  kcos = beta * cos(x - mu)
  constant = 2 * pi * b0
  if(log) d = kcos - log(constant) - beta
  else d = exp(kcos) / constant / exp(beta)
  d
}
dvm(1, 2, 20)
dvm(1:2, 2, 20)
dvm(1, 2:3, 20)
dvm(1:10, 1:2, 20)

dvm1(1, 2, 20, T)
dvm(1, 2, 20, T)
dvm1(1:10, 1, 20, T)
dvm(1:10, 1, 20, T)

dmixvm = function(x, mu=0, beta=1, pr=1, log=FALSE) {
  lm = length(mu)
  lx = length(x)
  if(length(beta) != lm) beta = rep(beta, length.out = lm)
  if(length(pr) != lm) pr = rep(pr, length.out = lm)  # beta&pr have same length as mu
  x.v = rep(x, lm)
  mu.v = rep(mu, rep(lx, lm))
  beta.v = rep(beta, rep(lx, lm))
  pr.v = rep(pr, rep(lx, lm))
  b0 = besselI(unique(beta.v), 0, expon.scaled = TRUE)
  b0.v = rep(rep(b0, length.out = lm), rep(lx, lm))
  kcos = beta.v * cos(x.v - mu.v)
  constant = 2 * pi * b0.v
  d = exp(kcos) / constant / exp(beta.v) * pr.v
  dim(d) = c(lx, lm)
  dmix = rowSums(d)
  if(log) dmix = log(dmix)
  dmix
}

dvmmix = function(x, mix, beta, log=FALSE) {
  n = length(x)
  m = length(mix$pt)
  logd = dvm(rep(x, m), rep(mix$pt, rep(n, m)), beta = beta, log = TRUE) + 
    rep(log(mix$pr), rep(n, m))
  dim(logd) = c(n, m)
  
  ma = apply(logd, 1, max)
  pid = exp(logd - ma)
  pis = rowSums(pid)
  r = ma + log(pis)
  if(!log) r = exp(r)
  r
}

(dvm(1:10, 1, 20) + dvm(1:10, 2, 20))/2
dmixvm(x = 1:10, mu = 1:2, beta = 20, pr = c(0.5, 0.5), FALSE)
dvmmix(1:10, disc(1:2, c(0.5, 0.5)), beta = 20, FALSE)
dmixvm(x = 1:10, mu = 1:2, beta = 20, pr = c(0.5, 0.5), TRUE)
dvmmix(1:10, disc(1:2, c(0.5, 0.5)), beta = 20, TRUE)

plot.npvm = function(x, mix, beta, nbins=50, dimen=1, area.prop=TRUE,
                     radius=1/sqrt(base::pi)) {
  m = length(mix$pt)
  z = seq(0, 2 * base::pi, len = 500)
  fd = function(x) dvmmix(x, mix, beta)
  d = fd(z)
  
  if(dimen == 1) {
    br = seq(0, 2 * base::pi, len = nbins + 1)
    whist(x$v, x$w, breaks = br, freq = FALSE, xlab = "Data", xaxt = "n",
          xlim = c(0, 2 * base::pi), ylim = range(0, d * 1.1),
          main=substitute("npvm (" * beta~"="~ a * "," ~ m~"="~b * ")",
                          list(a = signif(beta, 3), b = length(mix$pt))))
    axis(1, at = seq(0, 2, 0.5) * pi, padj = 0.5,
         labels = expression(0, frac(pi,2), pi, frac(3*pi,2), 2*pi))
    lines(z, d, col = "red", lwd = 1)
    points(mix$pt, rep(0, length(mix$pt)), col = 2)
    segments(mix$pt, rep(0, m), y1 = mix$pr * max(d), col = 2, lwd = 3)
  }
  else {
    chist(x$v, nbins = nbins, nlabels = 0, radius = radius, area.prop = area.prop)
    cdensity(fd, add = TRUE, radius = radius, area.prop = area.prop)
    points(radius * cos(mix$pt), radius * sin(mix$pt), col = 2)
    if(area.prop) h2 = sqrt(2 * mix$pr * max(d) + radius^2)
    else h2 = mix$pr * max(d) + radius
    segments(x0 = radius * cos(mix$pt), y0 = radius * sin(mix$pt),
             x1 = h2 * cos(mix$pt), y1 = h2 * sin(mix$pt), col = 2, lwd = 3)
  }
}

plot.npvm(as.npvm(seq(0, 2 * pi, len = 200)), disc(1:2, c(0.5, 0.5)), 20)

# mix = disc(pt = c(pi, pi/2), pr = c(0.3, 0.7))
# x = rnpvm(100, mu = mix$pt, pr = mix$pr, beta = 20, family = TRUE)
# plot(x, mix, beta = 20)
# cnm(x, init = list(beta = 10), plot = "p")

# ----



# Watson ------------------------------------------------------------------

rnpwa = function(n=100, mu=0, pr=1, beta=1, mix, family=FALSE) {
  if (missing(mix)) mix = disc(mu, pr)
  if (n == 0) return(numeric(0))
  l = length(mix$pt)
  if(length(beta) != l) beta = rep(beta, length.out = l)
  suppressWarnings(index <- sample.int(l, n, prob = mix$pr, replace = TRUE))
  mu.cart = matrix(c(cos(mu), sin(-mu)), nrow = 2, byrow = TRUE)   # polar to cartesian, ncol=no. components
  cart = list()
  polar = list()
  for (j in 1:l) {
    cart[[j]] = 
      Rfast::rbingham(sum(index == j), beta[j] * mu.cart[,j] %*% t(mu.cart[,j]))
    polar[[j]] = atan2(cart[[j]][,1], cart[[j]][,2])
    value = unlist(polar)
    value = ifelse(value < 0, value + pi, value)   # change the range to [0,pi]
  }
  x = list(v = value, w = 1)
  class(x) = "npwa"
  if (family) x
  else value
}

set.seed(17891642)
x1 = rnpwa(100, mu = c(pi, pi/2), pr = c(0.3, 0.7), beta = 20)

as.npwa = function(x, w=1) {
  x = list(v = c(x), w = w)
  class(x) = "npwa"
  x
}

npwa = function(v, w) {
  x = list(v = v, w = w)
  class(x) = "npwa"
  x
}

length.npwa = function(x) length(x$v)

weight.npwa = function(x, beta) x$w

gridpoints.npwa = function(x, beta, grid=100) {
  seq(0, base::pi, length = grid)
}

initial.npwa = function(x, beta, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = 20
  if(is.null(mix) || is.null(mix$pt)) {
    breaks = seq(0, base::pi, len = 20)
    r = nspmix::whist(x$v, x$w, breaks = breaks, freq = FALSE, plot = FALSE)
    i = r$density != 0
    mix = disc(r$mids[i], r$density[i])
  }
  list(beta = beta, mix = mix)
}

valid.npwa = function(x, beta, theta) beta > 0

suppspace.npwa = function(x, beta) c(0, base::pi)

logd.npwa = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld", "dk", "dt")

  if(which[1] == 1){
    dl$ld = beta * (cos(x$v- rep(pt, rep(length(x$v), length(pt))))^2) +
      log((base::pi * Re(fAsianOptions::kummerM(beta, 1/2, 1)))^(-1))
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  if(which[2] == 1)
    dl$db = NULL
  if(which[3] == 1){
    dl$dt = 2 * beta * sin(x$v - rep(pt, rep(length(x$v), length(pt)))) *
      cos(x$v - rep(pt, rep(length(x$v), length(pt))))
    dim(dl$dt) = c(length(x$v), length(pt))
  }
  dl
}

dwa1 = function(x=0, mu=0, beta=1, log=FALSE) {
  lx = length(x)
  lm = length(mu)
  x = rep(x, lm)
  mu = rep(mu, rep(lx, lm))
  m0 = Re(kummerM(beta, 1/2, 1, lnchf = 1))
  d = beta * (cos(x - mu))^2 - m0 - log(base::pi)
  dim(d) = c(lx, lm)
  if(log) d
  else exp(d)
}
dwa1(1, 2, 20)
dwa1(1:2, 2, 20)
dwa1(1, 2:3, 20)
dwa1(1:10, 1:2, 20)

dwa = function(x=0, mu=0, beta=1, log=FALSE) {
  lm = length(mu)
  lx = length(x)
  lb = length(beta)
  maxlen = max(lm, lx, lb)  # repeat to the largest length
  if(lx != maxlen) x = rep(x, length.out = maxlen)
  if(lm != maxlen) mu = rep(mu, length.out = maxlen)
  m0 = Re(fAsianOptions::kummerM(unique(beta), 1/2, 1, lnchf = 1))
  if(lb != maxlen) m0 = rep(m0, length.out = maxlen)
  d = beta * (cos(x - mu))^2 - m0 - log(base::pi)
  if(log) d
  else exp(d)
}
dwa(1, 2, 20)
dwa(1:2, 2, 20)
dwa(1, 2:3, 20)
dwa(1:10, 1:2, 20)

dwamix = function(x, mix, beta, log=FALSE) {
  n = length(x)
  m = length(mix$pt)
  logd = dwa(rep(x, m), rep(mix$pt, rep(n, m)), beta = beta, log = TRUE) + 
    rep(log(mix$pr), rep(n, m))
  dim(logd) = c(n, m)
  ma = apply(logd, 1, max)
  pid = exp(logd - ma)
  pis = rowSums(pid)
  r = ma + log(pis)
  if(!log) r = exp(r)
  r
}

(dwa(1:10, 1, 20) + dwa(1:10, 2, 20))/2
dwamix(1:10, disc(1:2, c(0.5, 0.5)), beta = 20, FALSE)
dwamix(1:10, disc(1:2, c(0.5, 0.5)), beta = 20, TRUE)

plot.npwa = function(x, mix, beta, nbins=50, dimen=1, area.prop=TRUE,
                     radius=1/sqrt(base::pi)) {
  m = length(mix$pt)
  z = seq(0, base::pi, len = 500)
  fd = function(x) dwamix(x, mix, beta)
  d = fd(z)
  
  if(dimen == 1) {
    br = seq(0, base::pi, len = nbins + 1)
    whist(x$v, x$w, breaks = br, freq = FALSE, xlab = "Data", xaxt = "n",
          xlim = c(0, base::pi), ylim = range(0, d * 1.1),
          main=substitute("npwa (" * beta~"="~ a * "," ~ m~"="~b * ")",
                          list(a = signif(beta, 3), b = length(mix$pt))))
    # axis(1, at = seq(0, 2, 0.5) * pi, padj = 0.5,
    #      labels = expression(0, frac(pi,2), pi, frac(3*pi,2), 2*pi))
    axis(1, at = seq(0, 1, 0.25) * pi, padj = 0.5,
         labels = expression(0, frac(pi,4), frac(pi,2), frac(3*pi,4), pi))
    lines(z, d, col = "red", lwd = 1)
    points(mix$pt, rep(0, length(mix$pt)), col = 2)
    segments(mix$pt, rep(0, m), y1 = mix$pr * max(d), col = 2, lwd = 3)
  }
  else {
    chist(x$v, nbins = nbins, nlabels = 0, radius = radius, area.prop = area.prop)
    cdensity(fd, add = TRUE, radius = radius, area.prop = area.prop)
    points(radius * cos(mix$pt), radius * sin(mix$pt), col = 2)
    if(area.prop) h2 = sqrt(2 * mix$pr * max(d) + radius^2)
    else h2 = mix$pr * max(d) + radius
    segments(x0 = radius * cos(mix$pt), y0 = radius * sin(mix$pt),
             x1 = h2 * cos(mix$pt), y1 = h2 * sin(mix$pt), col = 2, lwd = 3)
  }
}

plot.npwa(as.npwa(seq(0, pi, len = 200)), disc(1:2, c(0.5, 0.5)), 20)

# mix = disc(pt = c(pi, pi/2), pr = c(0.3, 0.7))
# x = rnpwa(100, mu = mix$pt, pr = mix$pr, beta = 20, family = TRUE)
# plot(x, mix, beta = 20)
# cnm(x, init = list(beta = 10), plot = "p")

# ----



# wrapped normal ----------------------------------------------------------

rnpwn = function(n=100, mu=0, pr=1, beta=0.5, mix, family=FALSE) {
  if (missing(mix)) mix = disc(mu, pr)
  if (n == 0) return(numeric(0))
  l = length(mix$pt)
  if(length(beta) != l) beta = rep(beta, length.out = l)
  suppressWarnings(index <- sample.int(l, n, prob = mix$pr, replace = TRUE))
  value = rep(0, n)
  for (i in 1:n) 
    value[i] = CircStats::rwrpnorm(1, mu = mix$pt[index[i]], rho = beta[index[i]])
  x = list(v=value, w=1)
  class(x) = "npwn"
  if (family) x  
  else value
}

set.seed(17891642)
x1 = rnpwn(100, mu = c(pi, pi/2), pr = c(0.3, 0.7), beta = 0.8)

as.npwn = function(x, w=1) {
  x = list(v=c(x), w=w)
  class(x) = "npwn"
  x
}

npwn = function(v, w) {
  x = list(v=v, w=w)
  class(x) = "npwn"
  x
}

length.npwn = function(x) length(x$v)

weight.npwn = function(x, beta) x$w

gridpoints.npwn = function(x, beta, grid=100) {
  seq(0, 2 * base::pi, length=grid)
}

gridpoints.npwn1 = function(x, beta, grid=100) {
  sd = sqrt(-2 * log(beta))
  breaks = pmax(ceiling(diff(range(x$v)) / (5*sd)), 5)   # number of breaks
  r = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
  i = r$density != 0
  i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
  m = sum(i)
  k = pmax(ceiling(grid / m), 10)
  d = r$breaks[2] - r$breaks[1]
  s = r$breaks[-length(r$breaks)][i]
  w = unique(rep(s, rep(k,m)) + d * 0:(k-1)/(k-1))
  w[w >=0 & w <= 2*pi]
}

initial.npwn = function(x, beta, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = 0.8
  if(is.null(mix) || is.null(mix$pt)) {
    sd = sqrt(-2 * log(beta))
    numbr = pmax(ceiling(diff(range(x$v)) / (5*sd)), 5)
    r = nspmix::whist(x$v, x$w, breaks=seq(0, 2*pi, length=numbr), freq=FALSE, plot=FALSE)
    i = r$density != 0
    mix = disc(r$mids[i], r$density[i])
  }
  list(beta=beta, mix=mix)
}

valid.npwn = function(x, beta, theta) beta >= 0 && beta <= 1

suppspace.npwn = function(x, beta) c(0, 2 * base::pi)

logd.npwn1 = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","dk","dt")

  # calculate logD
  var = -2 * log(beta)
  kmax = ceiling(sqrt(var)*10 / (2*pi))
  k = -kmax:kmax
  theta.l = rep(x$v, length(pt) * length(k))
  mu.l = rep(rep(pt, rep(length(x$v), length(pt))), length(k))
  k.l = rep(k, rep(length(x$v) * length(pt), length(k)))
  mat.s = theta.l - mu.l + 2 * k.l * pi
  dim(mat.s) = c(length(x$v) * length(pt), length(k))
  theta.s = rep(x$v, length(pt))
  mu.s = rep(pt, rep(length(x$v), length(pt)))
  d = theta.s - mu.s

  # use module to calculate the largest density
  d = ifelse(d < -pi, d + 4*pi, d)
  d = (d + pi) %% (2*pi) - pi
  d.large = d^2

  # log1p excluding the largest density
  d.small = exp((d.large - mat.s^2)/(2 * var))
  p = log1p(rowSums(d.small) - 1)
  dl$ld = -1/2 * log(2 * pi * var) - d.large / (2 * var) + p
  dim(dl$ld) = c(length(x$v), length(pt))

  if(which[1] == 1)
    dl$ld

  if(which[2] == 1)
    dl$db = NULL

  if(which[3] == 1){
    num = d.small * (mat.s / var)
    num.sum = rowSums(num)

    den = d.small
    den.sum = rowSums(den)

    dl$dt = num.sum / den.sum
    dim(dl$dt) = c(length(x$v), length(pt))
  }
  dl
}

# modified
logd.npwn2 = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","dk","dt")
  
  # calculate logD
  var = -2 * log(beta)
  kmax1 = max(1 / sqrt(var) * sqrt(-log(2 * pi^2 * var * er^2)), 
              sqrt(2) / sqrt(var))
  kmax2 = max(1 + sqrt(var) / pi * sqrt(-log(4 * pi^3 * er^2)), 
              1 / sqrt(var) * sqrt(-log(2 * pi^2 * var * er^2)))
  
  if(kmax2 > kmax1) {
    k = 1:ceiling(kmax2)
    theta.l = rep(x$v, length(pt) * length(k))
    mu.l = rep(rep(pt, rep(length(x$v), length(pt))), length(k))
    k.l = rep(k, rep(length(x$v) * length(pt), length(k)))
    
    term = cos(k.l * (theta.l - mu.l)) * beta^(k.l^2)
    dim(term) = c(length(x$v) * length(pt), length(k))
    dl$ld = log(2 * rowSums(term) + 1) - log(2 * pi)
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  
  else {
    kmax = ceiling(kmax1)
    k = -kmax:kmax
    theta.l = rep(x$v, length(pt) * length(k))
    mu.l = rep(rep(pt, rep(length(x$v), length(pt))), length(k))
    k.l = rep(k, rep(length(x$v) * length(pt), length(k)))
    
    mat.s = theta.l - mu.l + 2 * k.l * pi
    dim(mat.s) = c(length(x$v) * length(pt), length(k))
    theta.s = rep(x$v, length(pt))
    mu.s = rep(pt, rep(length(x$v), length(pt)))
    d = theta.s - mu.s
    
    # use module to calculate the largest density
    d = ifelse(d < -pi, d + 4*pi, d)
    d = (d + pi) %% (2*pi) - pi
    d.large = d^2
    
    # log1p excluding the largest density
    d.small = exp((d.large - mat.s^2)/(2 * var))
    p = log1p(rowSums(d.small) - 1)
    dl$ld = -1/2 * log(2 * pi * var) - d.large / (2 * var) + p
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  
  if(which[1] == 1)
    dl$ld
  
  if(which[2] == 1)
    dl$db = NULL
  
  if(which[3] == 1){
    if(kmax2 > kmax1) {
      num = 2 * k.l * beta^(k.l^2) * sin(k.l * (theta.l - mu.l))
      dim(num) = c(length(x$v) * length(pt), length(k))
      num.sum = rowSums(num)
      den.sum = 2 * rowSums(term) + 1
    }
    
    else{
      num = d.small * (mat.s / var)
      num.sum = rowSums(num)
      den = d.small
      den.sum = rowSums(den)
    }
    
    dl$dt = num.sum / den.sum
    dim(dl$dt) = c(length(x$v), length(pt))
  }
  dl
}

# modified again
logd.npwn = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","dk","dt")
  
  # calculate logD
  er = 1e-10
  var = -2 * log(beta)
  kmax1 = max(1 / sqrt(var) * sqrt(-log(2 * pi^2 * var * er^2)), 
              sqrt(2 / var))
  kmax2 = max(1 + sqrt(var) / pi * sqrt(-log(4 * pi^3 * er^2)), 
              1 + sqrt(var / 2) / pi)
  
  if(kmax2 > kmax1) {
    k = 1:ceiling(kmax2)
    theta.l = rep(x$v, length(pt) * length(k))
    mu.l = rep(rep(pt, rep(length(x$v), length(pt))), length(k))
    k.l = rep(k, rep(length(x$v) * length(pt), length(k)))
    
    term = cos(k.l * (theta.l - mu.l)) * beta^(k.l^2)
    dim(term) = c(length(x$v) * length(pt), length(k))
    dl$ld = log(2 * rowSums(term) + 1) - log(2 * pi)
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  
  else {
    kmax = ceiling(kmax1)
    k = -kmax:kmax
    theta.l = rep(x$v, length(pt) * length(k))
    mu.l = rep(rep(pt, rep(length(x$v), length(pt))), length(k))
    k.l = rep(k, rep(length(x$v) * length(pt), length(k)))
    
    term = theta.l - mu.l + 2 * k.l * pi
    dim(term) = c(length(x$v) * length(pt), length(k))
    dl$ld = log(rowSums(exp(-term^2 / 2 / var))) - log(2 * pi * var) / 2
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  
  if(which[1] == 1)
    dl$ld
  
  if(which[2] == 1)
    dl$db = NULL
  
  if(which[3] == 1){
    if(kmax2 > kmax1) {
      num = 2 * k.l * beta^(k.l^2) * sin(k.l * (theta.l - mu.l))
      dim(num) = c(length(x$v) * length(pt), length(k))
      num.sum = rowSums(num)
      den.sum = 2 * rowSums(term) + 1
    }                  
    
    else {
      num = exp(-term^2 / 2 / var) * term / var
      num.sum = rowSums(num)
      den = exp(-term^2 / 2 / var)
      den.sum = rowSums(den)
    }
    
    dl$dt = num.sum / den.sum
    dim(dl$dt) = c(length(x$v), length(pt))
  }
  dl
}

dwnmix = function(x, mix, beta) {
  d = suppressWarnings(outer(x, mix$pt, CircStats::dwrpnorm, rho = beta)) *
    rep(mix$pr, rep(length(x), length(mix$pt)))
  rowSums(d)
}

plot.npwn = function(x, mix, beta, nbins=50, dimen=1, area.prop=TRUE,
                     radius=1/sqrt(base::pi)) {
  m = length(mix$pt)
  z = seq(0, 2 * base::pi, len = 500)
  fd = function(x) dwnmix(x, mix, beta)
  d = fd(z)
  
  if(dimen == 1) {
    br = seq(0, 2 * base::pi, len = nbins + 1)
    whist(x$v, x$w, breaks = br, freq = FALSE, xlab = "Data", xaxt = "n",
          xlim = c(0, 2 * base::pi), ylim = range(0, d * 1.1),
          main=substitute("npwn (" * beta~"="~ a * "," ~ m~"="~b * ")",
                          list(a = signif(beta, 3), b = length(mix$pt))))
    axis(1, at = seq(0, 2, 0.5) * pi, padj = 0.5,
         labels = expression(0, frac(pi,2), pi, frac(3*pi,2), 2*pi))
    lines(z, d, col = "red", lwd = 1)
    points(mix$pt, rep(0, length(mix$pt)), col = 2)
    segments(mix$pt, rep(0, m), y1 = mix$pr * max(d), col = 2, lwd = 3)
  }
  else {
    chist(x$v, nbins = nbins, nlabels = 0, radius = radius, area.prop = area.prop)
    cdensity(fd, add = TRUE, radius = radius, area.prop = area.prop)
    points(radius * cos(mix$pt), radius * sin(mix$pt), col = 2)
    if(area.prop) h2 = sqrt(2 * mix$pr * max(d) + radius^2)
    else h2 = mix$pr * max(d) + radius
    segments(x0 = radius * cos(mix$pt), y0 = radius * sin(mix$pt),
             x1 = h2 * cos(mix$pt), y1 = h2 * sin(mix$pt), col = 2, lwd = 3)
  }
}

plot.npwn(as.npwn(x1), disc(1:2, c(0.5, 0.5)), 0.8)

mix = disc(pt = c(pi, pi/2), pr = c(0.3, 0.7))
x = rnpwn(100, mu = mix$pt, pr = mix$pr, beta = 0.8, family = TRUE)
plot(x, mix, beta = 0.8)
cnm(x, init = list(beta = 0.5), plot = "p")

# ----



# wrapped Cauchy ----------------------------------------------------------

rnpwc = function(n=100, mu=0, pr=1, beta=0.5, mix, family=FALSE) {
  if (missing(mix)) mix = disc(mu, pr)
  if (n == 0) return(numeric(0))
  l = length(mix$pt)
  if(length(beta) != l) beta = rep(beta, length.out = l)
  suppressWarnings(index <- sample.int(l, n, prob = mix$pr, replace = TRUE))
  value = rep(0, n)
  for (i in 1:n) 
    value[i] = CircStats::rwrpcauchy(1, location = mix$pt[index[i]], 
                                     rho = beta[index[i]])
  x = list(v=value, w=1)
  class(x) = "npwc"
  if (family) x  
  else value
}

set.seed(17891642)
x1 = rnpwc(100, mu = c(pi, pi/2), pr = c(0.3, 0.7), beta = 0.8)

as.npwc = function(x, w=1) {
  x = list(v=c(x), w=w)
  class(x) = "npwc"
  x
}

npwc = function(v, w) {
  x = list(v=v, w=w)
  class(x) = "npwc"
  x
}

length.npwc = function(x) length(x$v)

weight.npwc = function(x, beta) x$w

gridpoints.npwc = function(x, beta, grid=100) {
  seq(0, 2 * pi, length = grid)
}

initial.npwc = function(x, beta, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = 0.8
  if(is.null(mix) || is.null(mix$pt)) {
    breaks = seq(0, 2*pi, length=20)
    r = nspmix::whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
    i = r$density != 0
    mix = disc(r$mids[i], r$density[i])
  }
  list(beta=beta, mix=mix)
}

valid.npwc = function(x, beta, theta) beta >= 0 && beta <= 1

suppspace.npwc = function(x, beta) c(0, 2 * base::pi)

logd.npwc = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")

  if(which[1] == 1){
    dl$ld = log(CircStats::dwrpcauchy(theta = x$v, 
                                      mu = rep(pt,each=length(x$v)), rho=beta))
    dim(dl$ld) = c(length(x$v), length(pt))
  }
  if(which[2] == 1)
    dl$db = NULL
  if(which[3] == 1){
    dl$dt = 2*beta*sin(x$v-rep(pt,each=length(x$v))) / 
      (1+beta^2-2*beta*cos(x$v-rep(pt,each=length(x$v))))
    dim(dl$dt) = c(length(x$v), length(pt))
  }
  dl
}

dwcmix = function(x, mix, beta) {
  d = (outer(x, mix$pt, CircStats::dwrpcauchy, rho = beta)) *
    rep(mix$pr, rep(length(x), length(mix$pt)))
  rowSums(d)
}

plot.npwc = function(x, mix, beta, nbins=50, dimen=1, area.prop=TRUE,
                     radius=1/sqrt(base::pi)) {
  m = length(mix$pt)
  z = seq(0, 2 * base::pi, len = 500)
  fd = function(x) dwcmix(x, mix, beta)
  d = fd(z)
  
  if(dimen == 1) {
    br = seq(0, 2 * base::pi, len = nbins + 1)
    whist(x$v, x$w, breaks = br, freq = FALSE, xlab = "Data", xaxt = "n",
          xlim = c(0, 2 * base::pi), ylim = range(0, d * 1.1),
          main=substitute("npwc (" * beta~"="~ a * "," ~ m~"="~b * ")",
                          list(a = signif(beta, 3), b = length(mix$pt))))
    axis(1, at = seq(0, 2, 0.5) * pi, padj = 0.5,
         labels = expression(0, frac(pi,2), pi, frac(3*pi,2), 2*pi))
    lines(z, d, col = "red", lwd = 1)
    points(mix$pt, rep(0, length(mix$pt)), col = 2)
    segments(mix$pt, rep(0, m), y1 = mix$pr * max(d), col = 2, lwd = 3)
  }
  else {
    chist(x$v, nbins = nbins, nlabels = 0, radius = radius, area.prop = area.prop)
    cdensity(fd, add = TRUE, radius = radius, area.prop = area.prop)
    points(radius * cos(mix$pt), radius * sin(mix$pt), col = 2)
    if(area.prop) h2 = sqrt(2 * mix$pr * max(d) + radius^2)
    else h2 = mix$pr * max(d) + radius
    segments(x0 = radius * cos(mix$pt), y0 = radius * sin(mix$pt),
             x1 = h2 * cos(mix$pt), y1 = h2 * sin(mix$pt), col = 2, lwd = 3)
  }
}

plot.npwc(as.npwc(seq(0, pi, len = 200)), disc(1:2, c(0.5, 0.5)), 0.8)

mix = disc(pt = c(pi, pi/2), pr = c(0.3, 0.7))
x = rnpwc(100, mu = mix$pt, pr = mix$pr, beta = 0.8, family = TRUE)
plot(x, mix, beta = 0.8)
cnm(x, init = list(beta = 0.5), plot = "p")

# ----





