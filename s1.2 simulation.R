
mu1 = pi; pr1 = 1; beta1 = 20
x1 = rnpvm(n = 100, mu = mu1, pr = pr1, beta = beta1)

# MDE - modified cross validation -----------------------------------------

# split the data, return index of the test set
testind = function(n=500, k=10, random=TRUE) {
  if(k > n) stop("out of bounds")
  if(random) ind = sample(rep(1:k, length.out = n), size = n)  # random order
  else ind = sort(rep(1:k, length.out = n))
  split(1:n, ind)
}
testind(11, 3)
testind(11, 3, random = FALSE)
testind(5, 3)

# cv error function (continuous only)
errorfun = function(x, comp="npvm", tu=1, cvfold=10, cvrepe=10, cvseed=1789,
                    ictrain=FALSE, ictest=FALSE, nint=50,
                    cvloss=c("ISE","KL"), leaveone=FALSE) {
  n = length(x)
  error = 0
  if(leaveone) cvfold = n  # leave-one-out cv
  if(cvseed > 0) set.seed(cvseed)  # correlated sampling
  for (i in 1:cvrepe) {
    testlist = testind(n, cvfold)  # random order of data to split into train&test
    for (j in 1:cvfold) {
      testind = testlist[[j]]  # test set index
      test = as.npvm(x[testind])
      train = as.npvm(x[-testind])  # training set
      trainmodel = cnm(train, init = list(beta = tu))
      # calculate loss
      emp = dmixvm(test$v, trainmodel$mix$pt, trainmodel$beta, trainmodel$mix$pr)
      fsq = sqweight(trainmodel$mix, trainmodel$beta)
      cvloss = match.arg(cvloss)
      if(cvloss == "ISE") error = error + fsq - 2 * mean(emp)  # ISE
      else error = error - mean(log(emp))  # KL
    }
  }
  res = error / cvrepe / cvfold
  attributes(res) = list("numcomp" = length(trainmodel$mix$pt))
  res
}
errorfun(x1)
errorfun(x1, cvloss = "KL")

# cross validation
bws.cv = function(x, comp="npvm", cvseed=1789, inc=5e-03, numinc=5, opc=10,
                  lower=0.8, upper=1.2, seqlen=25, plotmodel=FALSE, plotpe=FALSE,
                  mlekappa=c("exact","approx"), hseq=NULL,
                  cvloss=c("ISE","KL","HD"), cvrepe=2, cvfold=10,
                  ictrain=FALSE, ictest=FALSE, nint=50,
                  print=FALSE, ploteach=FALSE) {
  # bandwidth sequence
  if(is.null(hseq)) {
    seq = c(seq(lower, 1, len = 3), upper^seq_len(pmax(1, seqlen - 3)))  # geometric seq, minimum len 4
    cirsd = switch(match.arg(mlekappa), "exact"=finvA(x), "approx"=mle.kappa(x))
    tu = cirsd * seq
  }
  else tu = hseq

  # set up initial values
  n = length(x)
  last = Inf
  allpe = rep(NA, length(tu))

  # for each tuning parameter
  for (i in 1:length(tu)) {
    pe = errorfun(x, tu = tu[i], cvseed = cvseed, cvloss = match.arg(cvloss),
                  cvrepe = cvrepe, cvfold = cvfold,
                  ictrain = ictrain, ictest = ictest, nint = nint)
    allpe[i] = pe
    d = diff(allpe) > 0  # whether current pe is larger than the preceeding one
    if(attributes(pe)$numcomp > n/opc) break  # stop if too many components
    # rle$length[rle$value] to extract only the TRUEs
    if (any(rle(d)$length[rle(d)$values] == numinc, na.rm = TRUE)) break  # stop if loss increase 3 times
    if(pe > last + inc) break  # stop if the increase in pe is larger than inc
    last = pe
    if(print) print(c(tu[i], pe, attributes(pe)$numcomp))
    if(ploteach) {  # plot for each bandwidth
      xvm = as.npvm(x)
      eachfit = cnm(xvm, init = list(beta = tu[i]), plot = "n")
      plot.npvm(xvm, eachfit$mix, eachfit$beta)
    }
  }
  # plot loss against bandwidth
  if(plotpe) {
    nonzero = !is.na(allpe)
    plot(tu[nonzero], allpe[nonzero], type = "o", pch = 16,
         xlab = "kappa", ylab = "loss", main = paste0("cvseed = ", cvseed))
  }

  # find bandwidth with the lowest pe
  betahat = tu[which.min(allpe)]
  res = cnm(as.npvm(x), init = list(beta = betahat), plot = "n")
  if(plotmodel) plot.npvm(as.npvm(x), res$mix, res$beta)
  res
}

mde.cv = bws.cv(x1, plotmodel = TRUE)
mde.cv
mde.cv.KL = bws.cv(x1, plotmodel = TRUE, cvloss = "KL")
mde.cv.KL
bws.cv(x1, plotmodel = TRUE, ictrain = TRUE)
bws.cv(x1, plotmodel = TRUE, ictest = TRUE)
bws.cv(x1, plotmodel = TRUE, ictrain = TRUE, ictest = TRUE, inc = 0.1)
bws.cv(x1, plotmodel = TRUE, ictrain = TRUE, ictest = TRUE, cvloss = "KL")

# ----



# MDE - AICc --------------------------------------------------------------

# based on mle.kappa
bws.AICc = function(x, comp="npvm", inc=3, lower=0.8, upper=1.2, seqlen=25,
                    plot="n", verb=0, plotmodel=FALSE, ctol=4, numinc=4, opc=10,
                    mlekappa=c("exact", "approx"), hseq=NULL,
                    censor=FALSE, nint=50, print=FALSE) {
  # bandwidth sequence
  if(is.null(hseq)) {
    seq = c(seq(lower, 1, len = 3), upper^seq_len(pmax(1, seqlen - 3)))  # geometric seq, minimum len 4
    cirsd = switch(match.arg(mlekappa), "exact"=finvA(x), "approx"=mle.kappa(x))
    tu = cirsd * seq
  }
  else tu = hseq

  if(is.null(nint)) nint = floor(length(data) + 10)  # number of intervals
  if(censor) x = fnpic(x, nint = nint)
  else x = as.npvm(x)

  # set up initial values
  n = length(x)
  last = Inf
  whole = list()
  allAICc = rep(NA, length(tu))

  # for each tuning parameter
  for (i in 1:length(tu)) {
    fit = cnm(x, init = list(beta = tu[i]), plot = plot, verbose = verb)
    fit$mix = nspmix:::collapse.snpmle(fit$mix, fit$beta, x, tol = ctol)
    whole[[i]] = fit
    loglh = fit$ll
    m = 2 * length(fit$mix$pt) - 1  # no. parameters
    AIC = 2 * m - 2 * loglh
    AICc = AIC + 2 * m * (m+1) / pmax(n - m - 1, 0)  # if no.para > no.data, den=0
    allAICc[i] = AICc
    d = diff(allAICc) > 0  # whether current AICc is larger than the preceeding one
    if (length(fit$mix$pt) > n/opc) break  # stop if too many components
    # rle$length[rle$value] to extract only the TRUEs
    if (any(rle(d)$length[rle(d)$values] == numinc, na.rm = TRUE)) break  # stop if AICc increase 3 times
    if (AICc > last + inc) break  # stop if AICc starts to increase by specified amount
    last = AICc
    if(print) print(c(tu[i], AICc, length(fit$mix$pt)))
  }

  # find the one with lowest AICc
  ind = which.min(allAICc)
  model = whole[[ind]]
  if(plotmodel) plot(x, model$mix, model$beta)
  # print(length(model$mix$pt))
  model
}

mde.AICc = bws.AICc(x1, plotmodel = TRUE)
mde.AICc
bws.AICc(x1, plotmodel = TRUE)
bws.AICc(x1, plotmodel = TRUE, censor = TRUE, nint = 50)

# ----



# KDE ---------------------------------------------------------------------

## circular package
?bw.cv.mse.circular  # squared error
?bw.cv.ml.circular   # Kullback-Leibler
?bw.nrd.circular     # rule-of-thumb

bw1 = bw.cv.mse.circular(circular(x1))
bw2 = bw.cv.ml.circular(circular(x1))
bw3 = bw.nrd.circular(circular(x1))
bw4 = bw.nrd.circular(circular(x1), kappa.est = "trigmoments")
c(bw1, bw2, bw3, bw4)

bw.cv.mse.circular(circular(x2), kernel = "wrappednormal")
bw.cv.ml.circular(circular(x2), kernel = "wrappednormal")
bw.nrd.circular(circular(x2))
bw.nrd.circular(circular(x2), kappa.est = "trigmoments")


## NPCirc package
?bw.boot
?bw.CV
?bw.pi
?bw.rt

bw5 = bw.boot(circular(x1))
bw6 = bw.CV(circular(x1))
bw7 = bw.pi(circular(x1))
bw8 = bw.rt(circular(x1))
c(bw5, bw6, bw7, bw8)

# ----



# tune the model and calculate loss  ----------------------------------------

# tune the model first
tune = function(method, bwcrt, data, comp="npvm", cvseed=1789, aiccinc=10,
                cvinc=5e-03, cvloss=c("ISE", "KL","HD"), cvrepe=2, cvfold=10,
                hseq=NULL, lower=0.8, upper=1.2, seqlen=25,
                mlekappa=c("exact","approx"), censor=FALSE,
                ictrain=FALSE, ictest=FALSE, nint=NULL) {
  if(method == "kde")
    switch(bwcrt,
           # "boot" = bw.boot(circular(data)),
           "LCV" = bw.CV(circular(data), method = "LCV"),
           "LSCV" = bw.CV(circular(data), method = "LSCV"),
           "pi" = bw.pi(circular(data)),
           "rt" = bw.rt(circular(data)))
  else
    switch(bwcrt,
           "AICc" = bws.AICc(data, comp = comp, mlekappa = mlekappa,
                             censor = censor, nint = nint, inc = aiccinc,
                             upper = upper, seqlen = seqlen, hseq = hseq),
           "cv" = bws.cv(data, comp = comp, cvseed = cvseed, inc = cvinc,
                         cvrepe = cvrepe, cvfold = cvfold,
                         upper = upper, seqlen = seqlen, cvloss = match.arg(cvloss),
                         mlekappa = match.arg(mlekappa), hseq = hseq,
                         ictrain=ictrain, ictest=ictest, nint=nint))
}

tune.pi = tune("kde", "pi", x1)
tune.AICc = tune("mde", "AICc", x1)
tune.cv = tune("mde", "cv", x1)
tune.cv.KL = tune("mde", "cv", x1, cvloss = "KL")
tune.AICc.ic = tune("mde", "AICc", x1, censor = TRUE)
tune.cv.ic = tune("mde", "cv", x1, ictest = TRUE)


# then calculate the loss
lossfun = function(loss=c("ISE","KL","HD"), method=c("kde","mde"), trainres,
                   data, simfam="npvm", mu=0, beta=1, alpha=0, pr=1) {
  n = length(data)

  dmix = switch(simfam, "npvm" = dmixvm, "npwn" = dmixwn, "npwc" = dmixwc,
                "npwsn" = dmixwsn)
  fden = function(x) dmix(x, mu, beta, pr)  # true density
  if(simfam == "npwsn") fden = function(x) dmix(x, mu, beta, alpha, pr)

  # if(simfam == "npvm") fden = function(x) dmixvm(x, mu, beta, pr)
  # else if(simfam == "npwn") fden = function(x) dmixwn(x, mu, beta, pr)
  # else if(simfam == "npwc") fden = function(x) dmixwc(x, mu, beta, pr)
  # else fden = function(x) dmixwsn(x, mu, beta, alpha, pr)

  fkde = function(x, bw) dmixvm(x, data, bw, rep(1 / n, n))  # est kde
  fmde = function(x, mix, beta) dmixvm(x, mix$pt, beta, mix$pr)  # est mde

  ISE.kde = function(x, bw) (fden(x) - fkde(x, bw))^2
  ISE.mde = function(x, mix, beta) (fden(x) - fmde(x, mix, beta))^2
  KL.kde = function(x, bw) fden(x) * log(fden(x) / fkde(x, bw))
  KL.mde = function(x, mix, beta) fden(x) * log(fden(x) / fmde(x, mix, beta))
  HD.kde = function(x, bw) (sqrt(fkde(x, bw)) - sqrt(fden(x)))^2
  HD.mde = function(x, mix, beta) (sqrt(fmde(x, mix, beta)) - sqrt(fden(x)))^2

  # ftrue2 = sqweight(disc(pt = mu, pr = pr), beta)
  # ftrue1 = integrate(function(x) fden(x)^2, 0, 2 * pi)$value
  # print(c(ftrue2, ftrue1))
  if(method == "kde") {
    switch(loss,
           # "ISE" = {
           #   fhat2 = sqweight(disc(pt = data, pr = rep(1/n, n)), trainres)
           #   fcross = crossweight(disc(pt = mu, pr = pr), beta,
           #                        disc(pt = data, pr = rep(1/n, n)), trainres)
           #   ftrue2 + fhat2 - 2 * fcross
           # },
           "ISE" = integrate(ISE.kde, 0, 2 * pi, bw = trainres)$value,
           "KL" = integrate(KL.kde, 0, 2 * pi, bw = trainres)$value,
           "HD" = integrate(HD.kde, 0, 2 * pi, bw = trainres)$value)
    # "ISE" = simpson1(ISE.kde, 0, 2 * pi, bw = trainres, n = 1000),
    # "KL" = simpson1(KL.kde, 0, 2 * pi, 100, bw = trainres, n = 1000),
    # "HD" = simpson1(HD.kde, 0, 2 * pi, bw = trainres, n = 1000)
  }
  else{
    switch(loss,
           # "ISE" = {
           #   fhat2 = sqweight(trainres$mix, trainres$beta)
           #   fcross = crossweight(disc(pt = mu, pr = pr), beta,
           #                        trainres$mix, trainres$beta)
           #   ftrue2 + fhat2 - 2 * fcross
           # },
           "ISE" = integrate(ISE.mde, 0, 2 * pi, mix = trainres$mix,
                             beta = trainres$beta)$value,
           "KL" = integrate(KL.mde, 0, 2 * pi, mix = trainres$mix,
                            beta = trainres$beta)$value,
           "HD" = integrate(HD.mde, 0, 2 * pi, mix = trainres$mix,
                            beta = trainres$beta)$value)
    # "ISE" = simpson1(ISE.mde, 0, 2 * pi, mix = trainres$mix,
    #                  beta = trainres$beta, n = 1000),
    # "KL" = simpson1(KL.mde, 0, 2 * pi, mix = trainres$mix,
    #                  beta = trainres$beta, n = 1000),
    # "HD" = simpson1(HD.mde, 0, 2 * pi, mix = trainres$mix,
    #                  beta = trainres$beta, n = 1000)
  }
}

lossfun(loss = "ISE", method = "kde", trainres = tune.pi, data = x1,
        mu = mu1, beta = beta1, pr = pr1)
lossfun(loss = "KL", method = "kde", trainres = tune.pi, data = x1,
        mu = mu1, beta = beta1, pr = pr1)
lossfun(loss = "ISE", method = "mde", trainres = tune.AICc, data = x1,
        mu = mu1, beta = beta1, pr = pr1)
lossfun(loss = "HD", method = "mde", trainres = tune.AICc, data = x1,
        mu = mu1, beta = beta1, pr = pr1)

# ----



# simulation function -----------------------------------------------------

### kde simulation
sim.kde = function(repetition=2, seed=1789, simfam="npvm",
                   n=200, mu=0, beta=1, alpha=0, pr=1, plot=TRUE) {
  # setup matrix to store results
  kdemethod = c("LCV", "LSCV", "pi", "rt")
  res = array(dim = c(repetition, 4, 3),
              dimnames = list(repe = 1:repetition, kdemethod = kdemethod,
                              loss = c("ISE", "KL", "HD")))

  # generate a set of random data for each repetition
  if (seed > 0) set.seed(seed)  # set seed to generate a series of random seed
  ranseed = sample.int(1e+08, repetition)
  for (i in 1:repetition) {
    set.seed(ranseed[i])
    if(simfam == "npvm") data = rnpvm(n=n, mu=mu, beta=beta, pr=pr)
    else if(simfam == "npwn") data = rnpwn(n=n, mu=mu, beta=beta, pr=pr)
    else if(simfam == "npwc") data = rnpwc(n=n, mu=mu, beta=beta, pr=pr)
    else data = rnpwsn(n=n, xi=mu, omega=beta, alpha=alpha, pr=pr)
    # print(head(data, 3))

    # 4 kde methods
    for(j in 1:4) {
      tuningbw = numeric(4)
      tuningbw[j] = tune(method="kde", bwcrt=j, data=data, comp="npvm")

      if(plot) {
        chist(data, nlabels = 0, nbins = 100,
              main = substitute(KDE ~ repetition ~ a ~ - b,
                                list(a = i, b = kdemethod[j])))
        dmix = switch(simfam, "npvm" = dmixvm, "npwn" = dmixwn, "npwc" = dmixwc,
                      "npwsn" = dmixwsn)
        ft = function(x) dmix(x, mu, beta, pr)
        if (simfam == "npwsn") ft = function(x) dmixwsn(x, mu, beta, alpha, pr)
        fe = function(x) dmixvm(x, data, tuningbw[j], rep(1/n, n))
        cdensity(ft, add = TRUE, col = 4)  # true density
        cdensity(fe, add = TRUE)           # estimated density
      }

      # 3 loss functions
      for (k in 1:3)
        res[i, j, k] = lossfun(loss=k, method="kde", trainres=tuningbw[j],
                               data=data, simfam=simfam,
                               mu=mu, beta=beta, alpha=alpha, pr=pr)
    }
  }

  list(mean = apply(res, 2:3, mean), se = apply(res, 2:3, sd) / repetition)
}

(qkde = sim.kde(2, mu = mu1, beta = beta1, pr = pr1, plot = TRUE))


### mde simulation with AICc
sim.AICc = function(repetition=2, seed=1789, simfam="npvm",
                    n=200, mu=0, beta=1, alpha=0, pr=1, aiccinc=100,
                    hseqfull=TRUE, hseq=NULL, lower=0.8, upper=1.2, seqlen=25,
                    censor=FALSE, nint=NULL, plot=TRUE) {
  # setup matrix to store results
  res = matrix(nrow = repetition, ncol = 3,
               dimnames = list(repe = 1:repetition, loss = c("ISE", "KL", "HD")))
  numcomp = numeric(repetition)
  betavalue = numcomp

  # generate a set of random data for each repetition
  if (seed > 0) set.seed(seed)  # set seed to generate a series of random seed
  ranseed = sample.int(1e+08, repetition)
  for (i in 1:repetition) {
    set.seed(ranseed[i])
    if(simfam == "npvm") data = rnpvm(n=n, mu=mu, beta=beta, pr=pr)
    else if(simfam == "npwn") data = rnpwn(n=n, mu=mu, beta=beta, pr=pr)
    else if(simfam == "npwc") data = rnpwc(n=n, mu=mu, beta=beta, pr=pr)
    else data = rnpwsn(n=n, xi=mu, omega=beta, alpha=alpha, pr=pr)
    # print(head(data, 3))

    # generate a bandwidth sequence based on full data
    if(hseqfull && is.null(hseq))
      hseq = fseq(finvA(data), lower=lower, upper=upper, seqlen=seqlen, plot=FALSE)

    trainmodel = tune(method="mde", bwcrt="AICc", data=data, comp="npvm",
                      hseq=hseq, lower=lower, upper=upper, seqlen=seqlen,
                      aiccinc=aiccinc, censor=censor, nint=nint)
    numcomp[i] = length(trainmodel$mix$pt)
    betavalue[i] = trainmodel$beta
    # print(length(trainmodel$mix$pt))

    if(plot) {
      chist(data, nlabels = 0, nbins = 100,
            main = substitute(AICc ~ repetition ~ a ~ b ~ c,
                              list(a = i,
                                   b = if(censor) paste0("bin", nint),
                                   c = paste0("m=", length(trainmodel$mix$pt)))))
      dmix = switch(simfam, "npvm" = dmixvm, "npwn" = dmixwn, "npwc" = dmixwc,
                    "npwsn" = dmixwsn)
      ft = function(x) dmix(x, mu, beta, pr)
      if (simfam == "npwsn") ft = function(x) dmixwsn(x, mu, beta, alpha, pr)
      fe = function(x) dmixvm(x, trainmodel$mix$pt,
                              trainmodel$beta, trainmodel$mix$pr)
      cdensity(ft, add = TRUE, col = 4)  # true density
      cdensity(fe, add = TRUE)           # estimated density
    }

    # 3 loss functions
    for (j in 1:3)
      res[i, j] = lossfun(loss=j, method="mde", trainres=trainmodel, data=data,
                          simfam=simfam, mu=mu, beta=beta, alpha=alpha, pr=pr)
  }
  list(mean = colMeans(res), se = apply(res, 2, sd) / repetition,
       numcomp = numcomp, betavalue = betavalue)
}

(qaicc = sim.AICc(2, mu = mu1, beta = beta1, pr = pr1, plot = TRUE))


### mde simulation with cv
sim.cv = function(repetition=2, seed=1789, simfam="npvm",
                  cvseed=1789, cvinc=1e-02, cvrepe=1, cvfold=10,
                  n=200, mu=0, beta=1, alpha=0, pr=1, plot=TRUE,
                  hseqfull=TRUE, hseq=NULL, lower=0.8, upper=1.2, seqlen=25,
                  mlekappa=c("exact","approx"), cvloss=c("ISE","KL"),
                  ictrain=FALSE, ictest=FALSE, nint=NULL) {
  # setup matrix to store results
  res = matrix(nrow = repetition, ncol = 3,
               dimnames = list(repe = 1:repetition, loss = c("ISE", "KL", "HD")))
  numcomp = numeric(repetition)
  betavalue = numcomp

  # generate a set of random data for each repetition
  if (seed > 0) set.seed(seed)  # set seed to generate a series of random seed
  ranseed = sample.int(1e+08, repetition)
  for (i in 1:repetition) {
    set.seed(ranseed[i])
    if(simfam == "npvm") data = rnpvm(n=n, mu=mu, beta=beta, pr=pr)
    else if(simfam == "npwn") data = rnpwn(n=n, mu=mu, beta=beta, pr=pr)
    else if(simfam == "npwc") data = rnpwc(n=n, mu=mu, beta=beta, pr=pr)
    else data = rnpwsn(n=n, xi=mu, omega=beta, alpha=alpha, pr=pr)
    # print(head(data, 3))

    # generate a bandwidth sequence based on full data
    if(hseqfull && is.null(hseq))
      hseq = fseq(finvA(data), lower=lower, upper=upper, seqlen=seqlen, plot=FALSE)

    trainmodel = tune(method="mde", bwcrt="cv", data=data, comp="npvm",
                      cvseed=cvseed, cvinc=cvinc, cvrepe=cvrepe,
                      cvloss=match.arg(cvloss),
                      hseq=hseq, lower=lower, upper=upper, seqlen=seqlen,
                      cvfold=cvfold, mlekappa=match.arg(mlekappa),
                      ictrain=ictrain, ictest=ictest, nint=nint)
    numcomp[i] = length(trainmodel$mix$pt)
    betavalue[i] = trainmodel$beta
    # print(length(trainmodel$mix$pt))

    if(plot) {
      chist(data, nlabels = 0, nbins = 100,
            main = substitute(CV ~ f ~ rep ~ a ~ b ~ c ~ d ~ e,
                              list(a = i,
                                   b = ifelse(ictrain, "ictrain", "") ,
                                   c = ifelse(ictest, "ictest", ""),
                                   d = if(ictrain | ictest == TRUE) paste0("bin", nint),
                                   e = paste0("m=", length(trainmodel$mix$pt)),
                                   f = match.arg(cvloss))))
      dmix = switch(simfam, "npvm" = dmixvm, "npwn" = dmixwn, "npwc" = dmixwc,
                    "npwsn" = dmixwsn)
      ft = function(x) dmix(x, mu, beta, pr)
      if (simfam == "npwsn") ft = function(x) dmixwsn(x, mu, beta, alpha, pr)
      fe = function(x) dmixvm(x, trainmodel$mix$pt,
                              trainmodel$beta, trainmodel$mix$pr)
      cdensity(ft, add = TRUE, col = 4)  # true density
      cdensity(fe, add = TRUE)           # estimated density
    }

    # 3 loss functions
    for (j in 1:3)
      res[i, j] = lossfun(loss = j, method = "mde", trainres = trainmodel,
                          data = data, simfam = simfam,
                          mu=mu, beta=beta, alpha=alpha, pr=pr)
  }

  list(mean = colMeans(res), se = apply(res, 2, sd) / repetition,
       numcomp = numcomp, betavalue = betavalue)
}

(qcv = sim.cv(2, mu = mu1, beta = beta1, pr = pr1, plot = TRUE, upper = 1.1))

# ----






