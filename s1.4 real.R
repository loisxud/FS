
# loss function for real data ---------------------------------------------

## tune does not need to change
tune("mde", "AICc", ant, censor = T)
tune("mde", "cv", ant, cvinc = 1)
tune("mde", "cv", ant, ictrain = T, ictest = T)

tunereal.lcv = tune("kde", "LCV", ant)
tunereal.AICc = tune("mde", "AICc", ant)
tunereal.cv = tune("mde", "cv", ant)
tunereal.AICc.ic = tune("mde", "AICc", ant, censor = TRUE)
tunereal.cv.ic = tune("mde", "cv", ant, ictest = TRUE)


## loss function for real data
# fullres
lossfun.real = function(fulldata, traindata, testdata, loss=c("ISE","KL"), 
                        method=c("kde","mde"), fullres, trainres) {
  nfull = length(fulldata)
  ntrain = length(traindata)
  if(method == "kde") {
    fn = dmixvm(testdata, traindata, trainres, rep(1/ntrain, ntrain))
    if(loss == "ISE") sqweight(disc(pt = fulldata, pr = rep(1/nfull, nfull)), 
                               fullres) - 2 * mean(fn)
    else -mean(log(fn))
  }
  else {
    fn = dmixvm(testdata, trainres$mix$pt, trainres$beta, trainres$mix$pr)
    if(loss == "ISE") sqweight(fullres$mix, fullres$beta) - 2 * mean(fn)
    else -mean(log(fn))
  }
}

# trainres only
lossfun.real = function(traindata, testdata, loss=c("ISE","KL"), 
                        method=c("kde","mde"), trainres) {
  ntrain = length(traindata)
  if(method == "kde") {
    fn = dmixvm(testdata, traindata, trainres, rep(1/ntrain, ntrain))
    if(loss == "ISE") sqweight(disc(pt = traindata, pr = rep(1/ntrain, ntrain)), 
                               trainres) - 2 * mean(fn)
    else -mean(log(fn))
  }
  else {
    fn = dmixvm(testdata, trainres$mix$pt, trainres$beta, trainres$mix$pr)
    if(loss == "ISE") sqweight(trainres$mix, trainres$beta) - 2 * mean(fn)
    else -mean(log(fn))
  }
}

lossfun.real(ant, ant, ant, "ISE", "mde", tunereal.AICc, tunereal.AICc)
lossfun.real(ant, ant, ant, "ISE", "kde", tunereal.lcv, tunereal.lcv)

lossfun.real(ant, ant, ant, "KL", "mde", tunereal.AICc, tunereal.AICc)
lossfun.real(ant, ant, ant, "KL", "kde", tunereal.lcv, tunereal.lcv)



## cross validation for real data
# mde
cv.real.mde = function(data, bwcrt=c("AICc","cv"), repetition=2, nfold=10, 
                       comp="npvm", seed=1789,  cvseed=1789, cvrepe=1, cvfold=10, 
                       plot=TRUE, aiccinc=100, cvinc=1e-02, cvloss=c("ISE","KL"),
                       hseqfull=TRUE, hseq=NULL, lower=0.8, upper=1.3, seqlen=25, 
                       censor=FALSE, ictrain=FALSE, ictest=FALSE, nint=NULL) {
  # setup matrix to store results
  res = array(dim = c(repetition, nfold, 2), 
              dimnames = list(repe = 1:repetition, nfold = 1:nfold,
                              loss = c("ISE", "KL")))
  numcomp = matrix(nrow = repetition, ncol = nfold, 
                   dimnames = list(repe = 1:repetition, nfold = 1:nfold))
  betavalue = numcomp
  n = length(data)
  bwcrt = match.arg(bwcrt)

  # generate a bandwidth sequence based on full data
  if(hseqfull && is.null(hseq))
    hseq = fseq(finvA(data), lower=lower, upper=upper, seqlen=seqlen, plot=FALSE)
  
  # tuning use full dataset first
  fullmodel = tune(method="mde", bwcrt=bwcrt, data=data, comp=comp, 
                   cvseed=cvseed, cvrepe=cvrepe, cvfold=cvfold, cvinc=cvinc, 
                   aiccinc=aiccinc, ictrain=ictrain, ictest=ictest,
                   hseq=hseq, lower=lower, upper=upper, seqlen=seqlen,
                   censor=censor, nint=nint, cvloss=match.arg(cvloss))
  
  # generate a set of random index for each repetition
  if(seed > 0) set.seed(seed)
  ranseed = sample.int(1e+08, repetition)
  for (i in 1:repetition) {
    # 10-fold cross validation
    set.seed(ranseed[i])
    testlist = testind(n, nfold)
    for (j in 1:nfold) {
      testind = testlist[[j]]  # test set index
      testset = data[testind]
      trainset = data[-testind]  # training set
      
      # construct model based on training data
      trainmodel = tune(method="mde", bwcrt=bwcrt, data=trainset, comp=comp, 
                        cvseed=cvseed, cvrepe=cvrepe, cvfold=cvfold, cvinc=cvinc, 
                        aiccinc=aiccinc, ictrain=ictrain, ictest=ictest,
                        hseq=hseq, lower=lower, upper=upper, seqlen=seqlen,
                        censor=censor, nint=nint, cvloss=match.arg(cvloss))
      numcomp[i, j] = length(trainmodel$mix$pt)
      betavalue[i, j] = trainmodel$beta
      
      if(plot) {
        chist(testset, nlabels = 0, nbins = 100, 
              main = substitute(e ~ f ~ repetition ~ a ~ fold ~ d ~ b ~ c, 
                                list(e = bwcrt, a = i, d = j,
                                     b = if(censor) paste0("bin", nint),
                                     c = paste0("m=", length(trainmodel$mix$pt)),
                                     f = if(bwcrt == "cv") cvloss)))
        fe = function(x) dmixvm(x, trainmodel$mix$pt, trainmodel$beta, 
                                trainmodel$mix$pr)
        cdensity(fe, add = TRUE)  # estimated density
      }
      
      # calculate 2 loss on test data
      for (k in 1:2)
        res[i, j, k] = lossfun.real(fulldata = data, traindata = trainset, testdata = testset,
                                    loss = c("ISE","KL")[k], method = "mde",
                                    fullres = fullmodel, trainres = trainmodel)
        # res[i, j, k] = lossfun.real(traindata = trainset, testdata = testset,
        #                             loss = c("ISE","KL")[k], method = "mde",
        #                             trainres = trainmodel)
    }
  }
  list(mean = c("ISE" = mean(res[, , 1]), "KL" = mean(res[, , 2])),
       se = c("ISE" = sd(res[, , 1]), "KL" = sd(res[, , 2])) / repetition / nfold, 
       numcomp = numcomp, betavalue = betavalue)
}

cv.real.mde(ant)
cv.real.mde(ant, bwcrt = "cv")
cv.real.mde(ant, bwcrt = "cv", cvloss = "KL")



# kde
cv.real.kde = function(data, repetition=2, nfold=10, comp="npvm", 
                       seed=1789, plot=TRUE) {
  # setup matrix to store results
  kdemethod = c("LCV", "LSCV", "pi", "rt")
  res = array(dim = c(repetition, nfold, 4, 2), 
              dimnames = list(repe = 1:repetition, nfold = 1:nfold,
                              kdemethod = kdemethod, loss = c("ISE", "KL")))
  n = length(data)
  
  # tuning use full dataset first
  fullbw = numeric(4)
  for (k in 1:4) fullbw[k] = tune(method="kde", bwcrt=k, data=data, comp="npvm")
  
  # generate a set of random index for each repetition
  if(seed > 0) set.seed(seed)
  ranseed = sample.int(1e+08, repetition)
  for (i in 1:repetition) {
    # 10-fold cross validation
    set.seed(ranseed[i])
    testlist = testind(n, nfold)
    for (j in 1:nfold) {
      testind = testlist[[j]]  # test set index
      testset = data[testind]
      trainset = data[-testind]  # training set
      
      # 4 kde methods
      for(k in 1:4) {
        trainbw = numeric(4)
        trainbw[k] = tune(method="kde", bwcrt=k, data=trainset, comp="npvm")
        
        if(plot) {
          chist(testset, nlabels = 0, nbins = 100, 
                main = substitute(KDE ~ repetition ~ a ~ fold ~ b ~ - c, 
                                  list(a = i, b = j, c = kdemethod[k])))
          fe = function(x) dmixvm(x, trainset, trainbw[k], 
                                  rep(1/length(trainset), length(trainset)))
          cdensity(fe, add = TRUE)  # estimated density
        }
        
        # 3 loss functions
        for (l in 1:2)
          res[i, j, k, l] = lossfun.real(fulldata=data, traindata=trainset, testdata=testset,
                                         loss=c("ISE","KL")[l], method="kde",
                                         fullres=fullbw[k], trainres=trainbw[k])
          # res[i, j, k, l] = lossfun.real(traindata=trainset, testdata=testset,
          #                                loss=c("ISE","KL")[l], method="kde",
          #                                trainres=trainbw[k])
      }
    }
  }    
  res.mean = apply(res, 3:4, mean)
  res.se = apply(res, 3:4, sd) / repetition / nfold
  list(mean = res.mean, se = res.se)
}

cv.real.kde(ant)

# ----



# real data ---------------------------------------------------------------

### CircNNTSR::Ants_radians (n = 100)
data(Ants_radians)
ant = Ants_radians
length(ant)
chist(ant, nlabels = 0, nbins = 100)

ant.kde = cv.real.kde(ant, repetition = 10)
ant.AICc = cv.real.mde(ant, "AICc", repetition = 10)
ant.cv = cv.real.mde(ant, "cv", repetition = 10)
ant.cv.KL = cv.real.mde(ant, "cv", repetition = 10, cvloss = "KL")

getres.real("ant")
getres.real("ant", "se")



### circular::fisherB1c (n = 254)
data(fisherB1c)
icu = unclass(fisherB1c)
attributes(icu) = NULL
icu = icu / 24 * 2 * pi
chist(icu, nlabels = 0, nbins = 100)

icu.kde = cv.real.kde(icu, repetition = 3)
icu.AICc = cv.real.mde(icu, "AICc", repetition = 10, seqlen = 20, upper = 1.5)
icu.cv = cv.real.mde(icu, "cv", repetition = 10, cvinc = 0.1, seqlen = 20, upper = 1.5)
icu.cv.KL = cv.real.mde(icu, "cv", repetition = 10, cvloss = "KL", cvinc = 0.1, seqlen = 20, upper = 1.5)

getres.real("icu")
getres.real("icu", "se")



### NPCirc::wind (n = 1752)
data(wind, package = "NPCirc")
length(wind$wind.dir)
winddir = wind$wind.dir
chist(winddir, nlabels = 0, nbins = 100)

winddir.kde = cv.real.kde(winddir, repetition = 3)
winddir.AICc = cv.real.mde(winddir, "AICc", repetition = 3, aiccinc = 100)
winddir.cv = cv.real.mde(winddir, "cv", repetition = 3)
winddir.cv.KL = cv.real.mde(winddir, "cv", repetition = 10, cvloss = "KL", 
                            cvinc = 0.1, seqlen = 20, upper = 1.5)

getres.real("winddir")
getres.real("winddir", "se")



# CircNNTSR::Turtles_radians
data(Turtles_radians)
turtle = Turtles_radians
length(turtle)
chist(turtle, nbins = 100, nlabels = 0)

f1 = function(x) dmixvm(x, turtle, 12.65736, rep(1/76,76))
cdensity(f1, add = TRUE)

a2 = bws.AICc(turtle, inc = 100)
f2 = function(x) dmixvm(x, a2$mix$pt, a2$beta, a2$mix$pr)
cdensity(f2, add = TRUE, col = 1)

turtle.kde = cv.real.kde(turtle, repetition = 3)
turtle.AICc = cv.real.mde(turtle, "AICc", repetition = 3)
turtle.cv = cv.real.mde(turtle, "cv", repetition = 3)
turtle.cv.KL = cv.real.mde(turtle, "cv", repetition = 10, cvinc = 0.1, seqlen = 20, upper = 1.5)

getres.real("turtle")
getres.real("turtle", "se")



# circular:fisherB6 (n = 100)
crossbed = unlist(fisherB6) / 180 * pi
length(crossbed)
chist(crossbed, nbins = 100, nlabels = 0)

crossbed.kde = cv.real.kde(crossbed, repetition = 3)
crossbed.AICc = cv.real.mde(crossbed, "AICc", repetition = 10, seqlen = 20, upper = 1.5)
crossbed.cv = cv.real.mde(crossbed, "cv", repetition = 10, cvinc = 0.1, 
                          seqlen = 20, upper = 1.5)
crossbed.cv.KL = cv.real.mde(crossbed, "cv", repetition = 10, cvloss = "KL", 
                             cvinc = 0.05, seqlen = 20, upper = 1.5)

getres.real("crossbed")
getres.real("crossbed", "se")
getres.real("crossbed", "numcomp")
getres.real("crossbed", "betavalue")

fseq(finvA(crossbed), upper = 1.25, seqlen = 30)
res = cnm(fnpvm(crossbed), init = list(beta = 5.3))
plot.npvm(fnpvm(crossbed), res$mix, res$beta, dimen = 2, nbins = 100)



# NPCirc::dragonfly (n = 214)
data(dragonfly)
dfly = dragonfly$orientation
length(dfly)
chist(dfly, nbins = 100, nlabels = 0)

dfly.kde = cv.real.kde(dfly, repetition = 3)
dfly.AICc = cv.real.mde(dfly, "AICc", repetition = 3, seqlen = 35, upper = 1.21)
dfly.cv = cv.real.mde(dfly, "cv", repetition = 3, cvinc = 0.1, seqlen = 35, upper = 1.21)
dfly.cv.KL = cv.real.mde(dfly, "cv", repetition = 3, cvloss = "KL", 
                         cvinc = 0.1, seqlen = 35, upper = 1.21)
dfly.cvone = cv.real.mde(dfly, "cv", repetition = 3, cvinc = 0.1, seqlen = 35, upper = 1.21)
dfly.cv.KLone = cv.real.mde(dfly, "cv", repetition = 3, cvloss = "KL", 
                            cvinc = 0.1, seqlen = 35, upper = 1.21)

getres.real("dfly")
getres.real("dfly", "se")
getres.real("dfly", "numcomp")
getres.real("dfly", "betavalue")

fseq(finvA(dfly), upper = 1.21, seqlen = 35)
res = cnm(fnpvm(dfly), init = list(beta = 105))
plot.npvm(fnpvm(dfly), res$mix, res$beta, dimen = 2, nbins = 100)



### CircNNTSR

# CircNNTSR hurricane
data(HurricanesGulfofMexico1951to1970)
hurricane = HurricanesGulfofMexico1951to1970$V1 * 2 * pi
length(hurricane)
chist(hurricane, nbins = 100)

hurricane.kde = cv.real.kde(hurricane, repetition = 3)
hurricane.AICc = cv.real.mde(hurricane, "AICc", repetition = 3)
hurricane.cv = cv.real.mde(hurricane, "cv", repetition = 3)

# CircNNTSR hurricane2
data(HurricanesGulfofMexico1971to2008)
hurricane2 = HurricanesGulfofMexico1971to2008$V1 * 2 * pi
length(hurricane2)
chist(hurricane2, nbins = 100)

hurricane2.kde = cv.real.kde(hurricane2, repetition = 3)
hurricane2.AICc = cv.real.mde(hurricane2, "AICc", repetition = 3, upper = 1.18, seqlen = 25)
hurricane2.cv = cv.real.mde(hurricane2, "cv", repetition = 3, upper = 1.18, seqlen = 25)
hurricane2.cv.KL = cv.real.mde(hurricane2, "cv", repetition = 3, cvloss="KL", 
                               upper = 1.18, seqlen = 25)

getres.real("hurricane2", mul = 10)

fseq(finvA(hurricane2), upper = 1.18, seqlen = 25)
res = cnm(fnpvm(hurricane2), init = list(beta = 5.6))
plot.npvm(fnpvm(hurricane2), res$mix, res$beta, dimen = 2, nbins = 100)



# CircNNTSR protein
data(ProteinsAAA)
dim(ProteinsAAA)
protein1 = ProteinsAAA$V1
protein2 = ProteinsAAA$V2
chist(protein1,  nbins = 100, nlabels = 0)
chist(protein2, nbins = 100, nlabels = 0)

fseq(finvA(protein1), upper = 1.14, seqlen = 30)
res = cnm(fnpvm(protein1), init = list(beta = 450))
length(res$mix$pt)
plot.npvm(fnpvm(protein1), res$mix, res$beta, dimen = 2, nbins = 100)

protein1.kde = cv.real.kde(protein1, repetition = 3)
protein1.AICc = cv.real.mde(protein1, "AICc", repetition = 3, upper = 1.14, seqlen = 30)
protein1.cv = cv.real.mde(protein1, "cv", repetition = 3, upper = 1.14, seqlen = 30)

getres.real("protein1")

protein2.kde = cv.real.kde(protein2, repetition = 3)
protein2.AICc = cv.real.mde(protein2, "AICc", repetition = 3)
protein2.cv = cv.real.mde(protein2, "cv", repetition = 3)

getres.real("protein2")



## large
# CircNNTSR wind
data(WindDirectionsTrivariate)
dim(WindDirectionsTrivariate)
chist(WindDirectionsTrivariate$V2, nbins = 100)
chist(WindDirectionsTrivariate$V6, nbins = 100)
chist(WindDirectionsTrivariate$V8, nbins = 100)



### NPCirc

# NPCirc::cross.beds1
data(cross.beds1)
beds = cross.beds1$angles
chist(beds, nbins = 100, nlabels = 0)

beds.kde = cv.real.kde(beds, repetition = 3)
beds.AICc = cv.real.mde(beds, "AICc", repetition = 3)
beds.cv = cv.real.mde(beds, "cv", repetition = 3)

rbind(AICc = beds.AICc$mean, beds.kde$mean)
rbind(AICc = beds.AICc$se, beds.kde$se)

# NPCirc::speed.wind2
data(speed.wind2)
speed2 = speed.wind2$Direction / 180 * pi
chist(speed2, nbins = 100, nlabels = 0, radius = 0.5)

# NPCirc::speed.wind
data(speed.wind)
dim(speed.wind)
speed = speed.wind$Direction / 180 * pi
speed = speed[!is.na(speed)]
chist(speed, nbins = 100, nlabels = 0, radius = 0.5)

# NPCirc::temp.wind
data(temp.wind)
dim(temp.wind)
temp = temp.wind$Direction / 180 * pi
temp = temp[!is.na(temp)]
chist(temp, nbins = 100, nlabels = 0)

temp.kde = cv.real.kde(temp, repetition = 1)
temp.AICc = cv.real.mde(temp, "AICc", repetition = 1)
temp.cv = cv.real.mde(temp, "cv", repetition = 1)

rbind(AICc = temp.AICc$mean, temp.kde$mean)
rbind(AICc = temp.AICc$se, temp.kde$se)

# ----



# simulated data ----------------------------------------------------------

# ordinary simulation
## 2 components
n2 = 100; mu2 = pi * c(0.5, 1.5); beta2 = 5; pr2 = c(1, 2)/3
chist(rnpvm(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixvm(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

# n = 100
xvm2.100.kde = sim.kde(3, n=100, mu=mu2, beta=beta2, pr=pr2)
xvm2.100.AICc = sim.AICc(3, n=100, mu=mu2, beta=beta2, pr=pr2)
xvm2.100.cvISE = sim.cv(3, n=100, mu=mu2, beta=beta2, pr=pr2)
xvm2.100.cvKL = sim.cv(3, n=100, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("xvm2", 100)
getres("xvm2", 100, "se")
getres("xvm2", 100, "betavalue")

# use sim data as real data
## 2 components
n2 = 100; mu2 = pi * c(0.5, 1.5); beta2 = 5; pr2 = c(1, 2)/3
xvm2 = rnpvm(n2, mu2, beta2, pr2)
chist(xvm2, nlabels = 0, nbins = 100)
f2 = function(x) dmixvm(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

fseq(finvA(xvm2), upper = 1.2, seqlen = 25)
res = cnm(fnpvm(xvm2), init = list(beta = 5.3))
plot.npvm(fnpvm(xvm2), res$mix, res$beta, dimen = 2, nbins = 100)

xvm2.kde = cv.real.kde(xvm2, repetition = 3)
xvm2.AICc = cv.real.mde(xvm2, "AICc", repetition = 3, upper = 1.2, seqlen = 25)
xvm2.cv = cv.real.mde(xvm2, "cv", repetition = 3, upper = 1.2, seqlen = 25)
xvm2.cv.KL = cv.real.mde(xvm2, "cv", repetition = 3, upper = 1.2, seqlen = 25,
                         cvloss = "KL")

getres.real("xvm2")
getres.real("xvm2", "se")
getres.real("xvm2", "betavalue")

## compare
summary(getres("xvm2", 100, "betavalue")[, 3])
summary(getres.real("xvm2", "betavalue")[, 3])

# ----





