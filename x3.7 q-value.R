
# sin^2 density -----------------------------------------------------------

# bm class
dbm.nor = function(x, mu=c(1,0,0), k=1, q=2, log=FALSE) {
  if(is.vector(x)) x = t(x)
  if(is.vector(mu)) mu = t(mu)
  p = ncol(x)
  if(ncol(mu) != p) stop("Dimension does not match.")
  x = x / sqrt(rowSums(x^2))
  mu = mu / sqrt(rowSums(mu^2))
  lx = nrow(x)
  lm = nrow(mu)
  lb = length(k)
  maxlen = max(lm, lx, lb)  # repeat to the largest length
  if(lx != maxlen) x = x[rep(1:lx, len = maxlen), , drop = FALSE]
  if(lm != maxlen) mu = mu[rep(1:lm, len = maxlen), , drop = FALSE]
  if(lb != maxlen) k = rep(k, length.out = maxlen)
  
  logd = -(k * (rowSums(mu * x)^2 - 1))^q
  logd = (-1)^(q+1) * (k * (rowSums(mu * x)^2 - 1))^q
  if(log) logd
  else exp(logd)
}

sort.bmmix = function(mix, by=1) {
  if(is.null(mix)) return(mix)
  index = order(mix$mu[, by])
  mix$mu = mix$mu[index, , drop = FALSE]
  mix$pr = mix$pr[index]
  if(length(mix$k) > 1) mix$k = mix$k[index]
  mix
}

bmmix = function(mu=c(1,0,0), pr=1, k=1, q=2, sort=TRUE, by=1) {
  if(is.vector(mu)) mu = t(mu)
  m = nrow(mu)
  p = ncol(mu)
  mu = mu / sqrt(rowSums(mu^2))
  # neg = mu[, p] < 0
  # mu[neg, ] = -mu[neg, ]
  if(length(pr) != m) pr = rep(pr, len = m)
  pr = pr / sum(pr)
  if(length(k) != m & length(k) > 1) pr = rep(pr, length.out = m)
  mix = list(mu = mu, pr = pr, k = k, q = q)
  class(mix) = "bmmix"
  if(sort) sort.bmmix(mix, by = by)
  else mix
}

outer.dbm.nor = function(x, mu, k, q=2, log=FALSE) {
  if(is.vector(x)) x = t(x)
  if(is.vector(mu)) mu = t(mu)
  n = nrow(x)
  m = nrow(mu)
  k = rep(k, length.out = m)
  d = dbm.nor(x = x[rep(1:n, m), , drop = FALSE], 
              mu = mu[rep(1:m, rep(n, m)), , drop = FALSE], 
              k = k[rep(1:m, rep(n, m))], q = q, log = log)
  dim(d) = c(n, m)
  d
}

dbmmix.nor = function(x, mix, log=FALSE) {
  if(is.vector(x)) x = t(x)
  j = mix$pr == 0
  if(any(j)) {
    mix$mu = mix$mu[!j, , drop = FALSE]
    mix$pr=mix$pr[!j]
  }
  logd = outer.dbm.nor(x, mix$mu, mix$k, mix$q, log = TRUE) + 
    rep(log(mix$pr), rep(nrow(x), length(mix$pr)))
  ma = apply(logd, 1, max)
  pid = exp(logd - ma)
  pis = rowSums(pid)
  r = ma + log(pis)
  if(!log) r = exp(r)
  r
}

logLik.bmmix.nor = function(mix, x, w=1) {
  d = dbmmix.nor(x, mix, log = TRUE)
  if(all(w == 1)) r = sum(d) 
  else r = sum(w * d)
  r
}

collapse.bmmix.nor = function(mix, x, coltol=1e-05, by=1) {
  mix = sort.bmmix(mix, by)
  j = mix$pr != 0  # remove components with zero mixing proportion
  mix$mu = mix$mu[j, , drop = FALSE]
  mix$pr = mix$pr[j]
  if(length(mix$pr) == 1) return(mix)
  ll = logLik.bmmix.nor(mix, x)
  mixt = mix
  
  for(i in 1:(length(mix$pr)-1)) {
    # compute pairwise cosine similarities of mixt$mu as distance measure
    # find index of the pair with the largest cosine similarity
    csmat = tcrossprod(mixt$mu)     # cosine similarities
    csmat[!upper.tri(csmat)] = NA   # remove duplicates and diagonal
    ip = c(arrayInd(which.max(abs(csmat)), dim(csmat)))  # index of the pair
    
    # remove the one with smaller mixing proportion
    irm = which.min(mixt$pr[ip])
    mixt$pr[ip[-irm]] = sum(mixt$pr[ip])
    mixt$pr = mixt$pr[-ip[irm]]
    mixt$mu = mixt$mu[-ip[irm], , drop = FALSE]
    
    llt = logLik.bmmix.nor(mixt, x)
    if(llt <= ll - coltol) break
    # if((llt-ll)/abs(ll) <= coltol) break
    else {ll = llt; mix = mixt}
  }
  # neg = mix$mu[, ncol(mix$mu)] < 0
  # mix$mu[neg, ] = -mix$mu[neg, ]
  # sort.bmmix(mix, by = by)
  mix
}

cnm.bm.nor = function(x, mix, k=1, q=2, maxlen=NULL, maxit=100, by=NULL, 
                      verbose=0, tol=1e-05, coltol=1e-05, pch=16, col=1, 
                      nlevels=5, extra=FALSE) {
  n = nrow(x)
  p = ncol(x)
  if(is.null(by)) by = which.max(apply(x, 2, sd))  # axis with the largest variation
  if(missing(mix) || is.null(mix)) {
    if(is.null(maxlen) || maxlen >= n) mix = bmmix(mu = x, pr = 1/n, k = k, q = q)
    else {
      ix = order(x[, by])
      sp = x[ix, , drop = FALSE]
      lx = nrow(x)
      ind = unique(round(seq(1, lx, len = maxlen)))
      mix = bmmix(mu = sp[ind, , drop = FALSE], pr = rep(1/lx, lx), k = k, q = q)
    }
  }
  else mix = sort.bmmix(mix, by = by)
  m = length(mix$pr)
  # neg = mix$mu[, p] < 0
  # mix$mu[neg, ] = -mix$mu[neg, ]
  ll = logLik.bmmix.nor(mix, x)
  convergence = 1
  if(maxit < 1)
    return(structure(list(num.iterations=0, convergence=1, max.gradient=NULL,
                          ll=ll, mix=mix)))
  
  order = order(mix$mu[, by])
  allsup = mix$mu[order, , drop = FALSE]  # sort in ascending order
  allden = outer.dbm.nor(x, allsup, mix$k, mix$q, log = TRUE)
  for (i in 1:maxit) {
    llt = ll
    if(verbose > 0) print.verbose(verbose, i-1, ll, mix)
    # neg = mix$mu[, p] < 0
    # mix$mu[neg, ] = -mix$mu[neg, ]
    
    # partition one axis at midpoints of consecutive support points
    mids = (mix$mu[-1, by] + mix$mu[-length(mix$pr), by]) * 0.5
    dmix = dbmmix.nor(x, mix, log = TRUE)
    g = colSums(exp(allden - dmix)) - n  # gradient
    
    # choose the support pt with the largest gradient in each interval
    index = indx(allsup[, by], mids) + 1
    jj = aggregate(g, by = list(group = index), which.max)  # position in each interval
    j = match(jj$group, index) + jj$x - 1  # position in the whole support set
    pt = allsup[j, , drop = FALSE]
    gpt = g[j]  # no need to check whether gradient is positive
    
    ## CNM
    mix2 = bmmix(mu = rbind(mix$mu, pt), 
                 pr = c(mix$pr, numeric(length(pt))), k=mix$k, q=mix$q, by=by)
    d = outer.dbm.nor(x, mix2$mu, mix2$k, mix2$q, log = TRUE)
    S = exp(d - dmix)
    grad = colSums(S)
    r = pnnls(S, 2, sum = 1)
    sol = r$x / sum(r$x)
    mix = mix2
    mix$pr = sol
    dea = sum(grad * (mix$pr - mix2$pr)) * 0.333
    mixt = mix
    alpha = 1
    for(j in 1:50) {  # Armijo line search
      ll1 = logLik.bmmix.nor(mixt, x)
      if(ll1 >= ll + dea * alpha) break
      alpha = alpha * 0.5
      mixt$pr = (1 - alpha) * mix2$pr + alpha * mix$pr
    }
    
    # mix = zero.bmmix(mixt, tol = coltol, by = by)
    mix = collapse.bmmix.nor(mixt, x, coltol = coltol, by = by)
    ll = logLik.bmmix.nor(mix, x)
    if(ll <= llt + tol) {convergence = 0; break}  # stopping criterion
  }
  
  if(extra) plot(z1[1:i], type = "l", main = "log-likelihood", 
                 xlab = "iteration", ylab = "loglh")
  if(verbose > 0) print.verbose(verbose, i, ll, mix)
  structure(list(num.iterations = i, convergence = convergence,
                 max.gradient = max(gpt), ll = ll, mix = mix),
            class="cnmmb")
}

fs.bm.nor = function(x, bwseq, q=2, start=2, lb=10, factor=2, tol=1e-05, 
                     coltol=1e-05, exhaustive=FALSE, by=NULL, ratio=1,
                     maxlen=NULL, maxit=100, verbose=0, 
                     pch=16, col=1, nlevels=5, lty=1, lwd=1) {
  if(missing(bwseq)) {
    seq = c(seq(0.8, 1, len = 3),  factor^seq_len(max(1, lb - 3)))
    bwseq = start * seq    
  }
  else bwseq = unique(sort(bwseq, decreasing = FALSE))
  lb = length(bwseq)
  n = nrow(x)
  p = ncol(x)
  info = matrix(nrow = lb, ncol = 5, 
                dimnames = list(1:lb, c("bw", "no.comp", "no.iter", "max.grad", "ll")))
  allmix = list()
  for(i in 1:lb) {
    res = cnm.bm.nor(x, k=bwseq[i], q=q, tol=tol, coltol=coltol, maxlen=maxlen, 
                     maxit=maxit, by=by)
    info[i, ] = c(bwseq[i], length(res$mix$pr), res$num.iterations, 
                  res$max.gradient, res$ll)
    allmix[[i]] = res
    if(verbose >= 1)
      cat("Iter", i, ":\tbw = ", sprintf("%.3f", bwseq[i]), "\t#comp = ", 
          length(res$mix$pr), "\n", sep="")
    if(verbose >= 2) print(res$mix)
    if(!exhaustive && length(res$mix$pr) > n/ratio) break
  }
  list(info = info[1:i, ], allmix = allmix)
}



set.seed(17891642)
x = rnpvmf(1000, mu = diag(10), k = 20)

res = cnm.bm.nor(x, k = 2, q = 1, tol = 0, verb = 1)
res$max.gradient

res = cnm.bm.nor(x, k = 2, q = 3, tol = 0, verb = 1)
res$max.gradient

res = cnm.bm.nor(x, k = 2, q = 6, tol = 0, verb = 1)
res$max.gradient

res = cnm.bm.nor(x, k = 2, q = 8, tol = 0, verb = 1)
res$max.gradient

res = cnm.bm.nor(x, k = 2, q = 4, tol = 0, coltol = 0, verb = 1)
res$max.gradient

# ----



# musk data setup ---------------------------------------------------------

data(musk)
# muskm <- ksvm(Class~.,data=musk,kernel="rbfdot",C=1000)
# muskmdim(musk)
dim(musk)
colnames(musk)
table(musk$Class)

musk.full = musk
musk = musk.full[, 1:166]
set.seed(17891642); testprop = 0.2
testind = c(sample.int(207, 207 * testprop), 
            sample.int(269, 269 * testprop) + 207)
musk.train = t(as.matrix(musk[-testind, ]))

# standardise features
musk.stand = stfea(musk.train)
dim(musk.stand)
# PCA of observations
musk.pca = prcomp(musk.stand, center = FALSE, scale. = FALSE)$x
dim(musk.pca)
sum(musk.pca[10, ]^2)

# prediction using all features
set.seed(17891642); testprop = 0.2
testind = c(sample.int(207, 207 * testprop), 
            sample.int(269, 269 * testprop) + 207)
train.all = musk.full[-testind, ]
test.all = musk.full[testind, ]
svm.all = svm(Class ~ ., data = train.all, kernel = "polynomial", cost = 10, scale = FALSE)
pred.svm.all = predict(svm.all, test.all)
allfea.musk = mean(test.all$Class == pred.svm.all)  # 0.946

# ----



# tuning q-value ----------------------------------------------------------

### q = 2
fseq(start = 2, factor = 1.12, lb = 30)
musk.mat.pca.bm2 = fs.bm.nor(musk.pca, start = 2, factor = 1.12, lb = 30, 
                             q = 2, tol = 1e-10, verb = 1)
musk.mat.pca.bm2$info

# max.grad
musk.grad.pca.bm2 = numeric(30)
for(i in 1:30) musk.grad.pca.bm2[i] = musk.mat.pca.bm2$allmix[[i]]$max.gradient
boxplot(musk.grad.pca.bm2)

musk.acc.pca.bm2 = numeric(30)
for(i in 1:30) {
  res = musk.mat.pca.bm2$allmix[[i]]$mix
  cossim = abs(tcrossprod(res$mu, musk.pca))
  index.pca = apply(cossim, 1, which.max)
  train.fs = train.all[, c(167, index.pca)]
  test.fs = test.all[, c(167, index.pca)]
  svm.fs = svm(Class ~ ., data = train.fs, kernel = "polynomial", cost = 10, scale = FALSE)
  pred.svm.fs = predict(svm.fs, test.fs)
  musk.acc.pca.bm2[i] = mean(test.fs$Class == pred.svm.fs)  # 0.936
}
plot(musk.mat.pca.bm2$info[, 2], musk.acc.pca.bm2, type = "o", xlab = "#comp")
abline(h = allfea.musk, col = 2, lty = 2)



### q = 4
fseq(start = 2, factor = 1.07, lb = 30)
musk.mat.pca.bm4 = fs.bm.nor(musk.pca, bwseq = fseq(start = 2, factor = 1.07, lb = 30), 
                              q = 4, tol = 1e-10, verb = 1, exhaustive = TRUE)
musk.mat.pca.bm4$info

# max.grad
musk.grad.pca.bm4 = numeric(30)
for(i in 1:30) musk.grad.pca.bm4[i] = musk.mat.pca.bm4$allmix[[i]]$max.gradient
boxplot(musk.grad.pca.bm4)

musk.acc.pca.bm4 = numeric(30)
for(i in 1:30) {
  res = musk.mat.pca.bm4$allmix[[i]]$mix
  cossim = abs(tcrossprod(res$mu, musk.pca))
  index.pca = apply(cossim, 1, which.max)
  train.fs = train.all[, c(167, index.pca)]
  test.fs = test.all[, c(167, index.pca)]
  svm.fs = svm(Class ~ ., data = train.fs, kernel = "polynomial", cost = 10, scale = FALSE)
  pred.svm.fs = predict(svm.fs, test.fs)
  musk.acc.pca.bm4[i] = mean(test.fs$Class == pred.svm.fs)  # 0.936
}
plot(musk.mat.pca.bm4$info[, 2], musk.acc.pca.bm4, type = "o", xlab = "#comp")
abline(h = allfea.musk, col = 2, lty = 2)



### q = 6
fseq(start = 2, factor = 1.07, lb = 30)
musk.mat.pca.bm6 = fs.bm.nor(musk.pca, bwseq = fseq(start = 2, factor = 1.07, lb = 30), 
                             q = 6, tol = 1e-10, verb = 1, exhaustive = TRUE)
musk.mat.pca.bm6$info

# max.grad
musk.grad.pca.bm6 = numeric(30)
for(i in 1:30) musk.grad.pca.bm6[i] = musk.mat.pca.bm6$allmix[[i]]$max.gradient
boxplot(musk.grad.pca.bm6)

musk.acc.pca.bm6 = numeric(30)
for(i in 1:30) {
  res = musk.mat.pca.bm6$allmix[[i]]$mix
  cossim = abs(tcrossprod(res$mu, musk.pca))
  index.pca = apply(cossim, 1, which.max)
  train.fs = train.all[, c(167, index.pca)]
  test.fs = test.all[, c(167, index.pca)]
  svm.fs = svm(Class ~ ., data = train.fs, kernel = "polynomial", cost = 10, scale = FALSE)
  pred.svm.fs = predict(svm.fs, test.fs)
  musk.acc.pca.bm6[i] = mean(test.fs$Class == pred.svm.fs)  # 0.936
}
plot(musk.mat.pca.bm6$info[, 2], musk.acc.pca.bm6, type = "o", xlab = "#comp")
abline(h = allfea.musk, col = 2, lty = 2)



### q = 8
fseq(start = 2, factor = 1.07, lb = 30)
musk.mat.pca.bm8 = fs.bm.nor(musk.pca, bwseq = fseq(start = 2, factor = 1.07, lb = 30), 
                              q = 8, tol = 1e-10, verb = 1, exhaustive = TRUE)
musk.mat.pca.bm8$info

# max.grad
musk.grad.pca.bm8 = numeric(30)
for(i in 1:30) musk.grad.pca.bm8[i] = musk.mat.pca.bm8$allmix[[i]]$max.gradient
boxplot(musk.grad.pca.bm8)

musk.acc.pca.bm8 = numeric(30)
for(i in 1:30) {
  res = musk.mat.pca.bm8$allmix[[i]]$mix
  cossim = abs(tcrossprod(res$mu, musk.pca))
  index.pca = apply(cossim, 1, which.max)
  train.fs = train.all[, c(167, index.pca)]
  test.fs = test.all[, c(167, index.pca)]
  svm.fs = svm(Class ~ ., data = train.fs, kernel = "polynomial", cost = 10, scale = FALSE)
  pred.svm.fs = predict(svm.fs, test.fs)
  musk.acc.pca.bm8[i] = mean(test.fs$Class == pred.svm.fs)  # 0.936
}
plot(musk.mat.pca.bm8$info[, 2], musk.acc.pca.bm8, type = "o", xlab = "#comp")
abline(h = allfea.musk, col = 2, lty = 2)



### q = 10
fseq(start = 2, factor = 1.07, lb = 30)
musk.mat.pca.bm10 = fs.bm.nor(musk.pca, start = 2, factor = 1.07, lb = 30, 
                              q = 10, tol = 1e-10, verb = 1)
musk.mat.pca.bm10$info

# max.grad
musk.grad.pca.bm10 = numeric(30)
for(i in 1:30) musk.grad.pca.bm10[i] = musk.mat.pca.bm10$allmix[[i]]$max.gradient
boxplot(musk.grad.pca.bm10)

musk.acc.pca.bm10 = numeric(30)
for(i in 1:30) {
  res = musk.mat.pca.bm10$allmix[[i]]$mix
  cossim = abs(tcrossprod(res$mu, musk.pca))
  index.pca = apply(cossim, 1, which.max)
  train.fs = train.all[, c(167, index.pca)]
  test.fs = test.all[, c(167, index.pca)]
  svm.fs = svm(Class ~ ., data = train.fs, kernel = "polynomial", cost = 10, scale = FALSE)
  pred.svm.fs = predict(svm.fs, test.fs)
  musk.acc.pca.bm10[i] = mean(test.fs$Class == pred.svm.fs)  # 0.936
}
plot(musk.mat.pca.bm10$info[, 2], musk.acc.pca.bm10, type = "o", xlab = "#comp")
abline(h = allfea.musk, col = 2, lty = 2)



### q = 12
fseq(start = 2, factor = 1.07, lb = 30)
musk.mat.pca.bm12 = fs.bm.nor(musk.pca, bwseq = fseq(start = 2, factor = 1.07, lb = 30), 
                              q = 12, tol = 1e-10, verb = 1, exhaustive = TRUE)
musk.mat.pca.bm12$info

# max.grad
musk.grad.pca.bm12 = numeric(30)
for(i in 1:30) musk.grad.pca.bm12[i] = musk.mat.pca.bm12$allmix[[i]]$max.gradient
boxplot(musk.grad.pca.bm12)

musk.acc.pca.bm12 = numeric(30)
for(i in 1:30) {
  res = musk.mat.pca.bm12$allmix[[i]]$mix
  cossim = abs(tcrossprod(res$mu, musk.pca))
  index.pca = apply(cossim, 1, which.max)
  train.fs = train.all[, c(167, index.pca)]
  test.fs = test.all[, c(167, index.pca)]
  svm.fs = svm(Class ~ ., data = train.fs, kernel = "polynomial", cost = 10, scale = FALSE)
  pred.svm.fs = predict(svm.fs, test.fs)
  musk.acc.pca.bm12[i] = mean(test.fs$Class == pred.svm.fs)  # 0.936
}
plot(musk.mat.pca.bm12$info[, 2], musk.acc.pca.bm12, type = "o", xlab = "#comp")
abline(h = allfea.musk, col = 2, lty = 2)

# ----

### compare
plot(musk.mat.pca.bm2$info[, 2], musk.acc.pca.bm2, type = "o", xlab = "#comp", 
     ylim = c(0.5, 1))
abline(h = allfea.musk, col = 2, lty = 2)
lines(musk.mat.pca.bm4$info[, 2], musk.acc.pca.bm4, type = "o", col = 3, 
      ylim = c(0.5, 1))
lines(musk.mat.pca.bm6$info[, 2], musk.acc.pca.bm6, type = "o", col = 7, 
      ylim = c(0.5, 1))
lines(musk.mat.pca.bm8$info[, 2], musk.acc.pca.bm8, type = "o", col = 4, 
      ylim = c(0.5, 1))
lines(musk.mat.pca.bm10$info[, 2], musk.acc.pca.bm10, type = "o", col = 6, 
      ylim = c(0.5, 1))
lines(musk.mat.pca.bm12$info[, 2], musk.acc.pca.bm12, type = "o", col = 8, 
      ylim = c(0.5, 1))
lines(musk.mat.pca$info[, 2], musk.acc.pca, type = "o", col = "orange", 
      ylim = c(0.5, 1))  # Watson








