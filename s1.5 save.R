
savepath = "H:/Documents/2019/Visualisation/R.WIP/simulation/"
savepath2 = "H:/Documents/2019/Visualisation/R.WIP/realdata/"
savepath3 = "H:/Documents/2019/Visualisation/R.WIP/test/"

fbw = function(x, seqlen=25, upper=1.2, lower=0.8) {
  s = fseq(finvA(x), seqlen = seqlen, upper = upper, lower = lower)
  res = cnm(fnpvm(x), init = list(beta = s[length(s)]))
  plot.npvm(fnpvm(x), res$mix, res$beta, dimen = 2, nbins = 100)
  list("bandwidth sequence" = s, "number of component" = length(res$mix$pt))
}

# functions to aggregate results
getres = function(family="vm1", size=100, type=c("mean","se","numcomp","betavalue"), 
                  mul=100, digit=2) {
  AICc = get0(paste(family, size, "AICc", sep = "."))
  cvISE = get0(paste(family, size, "cvISE", sep = "."))
  cvKL = get0(paste(family, size, "cvKL", sep = "."))
  kde = get0(paste(family, size, "kde", sep = "."))
  
  type = match.arg(type)
  if(type == "mean") 
    round(rbind(AICc = AICc$mean, cvISE = cvISE$mean, cvKL = cvKL$mean, 
                kde$mean) * mul, digits = digit)
  else if(type == "se") 
    round(rbind(AICc = AICc$se, cvISE = cvISE$se, 
                cvKL = cvKL$se, kde$se) * mul, digits = digit)
  else if(type == "numcomp") 
    cbind(AICc = c(AICc$numcomp), cvISE = c(cvISE$numcomp), cvKL = c(cvKL$numcomp))
  else 
    round(cbind(AICc = c(AICc$betavalue), cvISE = c(cvISE$betavalue), 
                cvKL = c(cvKL$betavalue)), 2)
}
getres("vm1", 100)

getres.real = function(dataset="dfly", type=c("mean","se","numcomp","betavalue"),
                       mul=100, digit=2) {
  AICc = get0(paste(dataset, "AICc", sep = "."))
  cvISE = get0(paste(dataset, "cv", sep = "."))
  cvKL = get0(paste(dataset, "cv.KL", sep = "."))
  kde = get0(paste(dataset, "kde", sep = "."))
  
  type = match.arg(type)
  if(type == "mean") 
    round(rbind(AICc = AICc$mean, cvISE = cvISE$mean, cvKL = cvKL$mean, 
                kde$mean) * mul, digits = digit)
  else if(type == "se") 
    round(rbind(AICc = AICc$se, cvISE = cvISE$se, 
                cvKL = cvKL$se, kde$se) * mul, digits = digit)
  else if(type == "numcomp") 
    cbind(AICc = c(AICc$numcomp), cvISE = c(cvISE$numcomp), cvKL = c(cvKL$numcomp))
  else 
    round(cbind(AICc = c(AICc$betavalue), cvISE = c(cvISE$betavalue), 
                cvKL = c(cvKL$betavalue)), 2)
}
getres.real("turtle")

aggres = function(family="vm1", type=c("mean","se","numcomp","betavalue"), 
                  mul=100, digit=2) {
  type = match.arg(type)
  
  m100 = getres(family=family, size=100, type=type, mul=mul, digit=digit)
  m400 = getres(family=family, size=400, type=type, mul=mul, digit=digit)
  m1600 = getres(family=family, size=1600, type=type, mul=mul, digit=digit)
  m100 = m100[c(7:4, 1:3), ]
  m400 = m400[c(7:4, 1:3), ]
  m1600 = m1600[c(7:4, 1:3), ]
  
  array(c(m100, m400, m1600),
        dim = c(7, 3, 3), 
        dimnames = list("method" = c("RT","PI","LSCV","LCV","AICc","cvISE","cvKL"), 
                        "loss" = c("ISE", "KL", "HD"),
                        "size" = c(100, 400, 1600)))
}
aggres("vm1")

# convert to table in latex
fmin = function(x) which(x == min(x))
latexres = function(family="vm1", loss=c("ISE","KL","HD"), mul=100, digit=2) {
  loss = match.arg(loss)
  nloss = switch(loss, "ISE" = 1, "KL" = 2, "HD" = 3)
  
  mat = NULL
  for (i in 1:length(family)) {
    mean = aggres(family = family[i], type = "mean", mul = mul, digit = digit)[, nloss, ]
    se = aggres(family = family[i], type = "se", mul = mul, digit = digit)[, nloss, ]
    new = paste(format(c(mean), digits = digit),
                paste0("(", format(se, digits = digit), ")"), sep = " ")
    # ind = apply(mean, 2, which.min) + 7 * 0:2  # only one minimum
    ind = unlist(mapply("+", apply(mean, 2, fmin), 7 * 0:2))  # find multiple minimum
    new[ind] = paste0("\\textbf{", new[ind], "}")
    mat = cbind(mat, new)
  }
  c1 = c("KDE", rep("", 3), "MDE", rep("", 2))
  c2 = c("RT", "PI", "LSCV", "LCV", "AICc", "CVISE", "CVKL")
  mat = cbind(c1, c2, mat)
  result = apply(mat, 1, paste, collapse = " & ")
  result = paste0(result, "\\\\")
  res = character(length(result) + 10)
  res[-c(1:2, 7, 11:12, 17, 21:22, 27, 31)] = result
  res[1] = "\\toprule \\textit{Estimator} & \\textit{Selection method} & vM2 & WN3 & WC2 & WSN \\\\ \\toprule"
  res[c(2, 12, 22)] = paste0("&&\\multicolumn{4}{c}{$n = ", 
                             c(100,400,1600), 
                             "$} \\\\ \\cmidrule{3-6}")
  res[c(7, 17, 27)] = "\\midrule"
  res[c(11, 21, 31)] = "\\bottomrule"
  cat(res, fill = TRUE)
}

latexres(c("vm2", "wn3", "wc2", "wsn1"))
latexres(c("vm2", "wn3", "wc2", "wsn1"), loss = "KL")
latexres(c("vm2", "wn3", "wc2", "wsn1"), loss = "HD")



# von Mises mixture -------------------------------------------------------

## 1 component
n1 = 400; mu1 = pi * 0; beta1 = 20; pr1 = 1
chist(rnpvm(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixvm(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

x = rnpvm(100, mu1, beta1, pr1)
fbw(x, 40, 1.05)

# n = 100
pdf(file = paste0(savepath, "vm1.100.kde", ".pdf"))
vm1.100.kde = sim.kde(10, n=100, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "vm1.100.AICc", ".pdf"))
vm1.100.AICc = sim.AICc(10, n=100, mu=mu1, beta=beta1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath, "vm1.100.cvISE", ".pdf"))
vm1.100.cvISE = sim.cv(10, n=100, mu=mu1, beta=beta1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath, "vm1.100.cvKL", ".pdf"))
vm1.100.cvKL = sim.cv(10, n=100, mu=mu1, beta=beta1, cvloss="KL", seqlen=40, upper=1.05)
dev.off()

# n = 400
pdf(file = paste0(savepath, "vm1.400.kde", ".pdf"))
vm1.400.kde = sim.kde(10, n=400, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "vm1.400.AICc", ".pdf"))
vm1.400.AICc = sim.AICc(10, n=400, mu=mu1, beta=beta1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath, "vm1.400.cvISE", ".pdf"))
vm1.400.cvISE = sim.cv(10, n=400, mu=mu1, beta=beta1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath, "vm1.400.cvKL", ".pdf"))
vm1.400.cvKL = sim.cv(10, n=400, mu=mu1, beta=beta1, cvloss="KL", seqlen=40, upper=1.05)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "vm1.1600.kde", ".pdf"))
vm1.1600.kde = sim.kde(10, n=1600, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "vm1.1600.AICc", ".pdf"))
vm1.1600.AICc = sim.AICc(10, n=1600, mu=mu1, beta=beta1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath, "vm1.1600.cvISE", ".pdf"))
vm1.1600.cvISE = sim.cv(10, n=1600, mu=mu1, beta=beta1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath, "vm1.1600.cvKL", ".pdf"))
vm1.1600.cvKL = sim.cv(10, n=1600, mu=mu1, beta=beta1, cvloss="KL", seqlen=40, upper=1.05)
dev.off()



## 2 components
n2 = 400; mu2 = pi * c(0.2, 0.5); beta2 = 15; pr2 = c(1, 2)/3
chist(rnpvm(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixvm(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

x = rnpvm(100, mu2, beta2, pr2)
fbw(x, 40, 1.09)

# n = 100
pdf(file = paste0(savepath, "vm2.100.kde", ".pdf"))
vm2.100.kde = sim.kde(10, n=100, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "vm2.100.AICc", ".pdf"))
vm2.100.AICc = sim.AICc(10, n=100, mu=mu2, beta=beta2, pr=pr2, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm2.100.cvISE", ".pdf"))
vm2.100.cvISE = sim.cv(10, n=100, mu=mu2, beta=beta2, pr=pr2, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm2.100.cvKL", ".pdf"))
vm2.100.cvKL = sim.cv(10, n=100, mu=mu2, beta=beta2, pr=pr2, cvloss="KL", seqlen=40, upper=1.09)
dev.off()

# n = 400
pdf(file = paste0(savepath, "vm2.400.kde", ".pdf"))
vm2.400.kde = sim.kde(10, n=400, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "vm2.400.AICc", ".pdf"))
vm2.400.AICc = sim.AICc(10, n=400, mu=mu2, beta=beta2, pr=pr2, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm2.400.cvISE", ".pdf"))
vm2.400.cvISE = sim.cv(10, n=400, mu=mu2, beta=beta2, pr=pr2, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm2.400.cvKL", ".pdf"))
vm2.400.cvKL = sim.cv(10, n=400, mu=mu2, beta=beta2, pr=pr2, cvloss="KL", seqlen=40, upper=1.09)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "vm2.1600.kde", ".pdf"))
vm2.1600.kde = sim.kde(10, n=1600, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "vm2.1600.AICc", ".pdf"))
vm2.1600.AICc = sim.AICc(10, n=1600, mu=mu2, beta=beta2, pr=pr2, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm2.1600.cvISE", ".pdf"))
vm2.1600.cvISE = sim.cv(10, n=1600, mu=mu2, beta=beta2, pr=pr2, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm2.1600.cvKL", ".pdf"))
vm2.1600.cvKL = sim.cv(10, n=1600, mu=mu2, beta=beta2, pr=pr2, cvloss="KL", seqlen=40, upper=1.09)
dev.off()



### 3 components
n3 = 400; mu3 = pi * c(0.2, 0.5, 0.8); beta3 = 15; pr3 = c(1,3,2)/6
chist(rnpvm(n3, mu3, beta3, pr3), nlabels = 0, nbins = 100)
f3 = function(x) dmixvm(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE, col = 4)

x = rnpvm(100, mu3, beta3, pr3)
fbw(x, 40, 1.09)

# n = 100
pdf(file = paste0(savepath, "vm3.100.kde", ".pdf"))
vm3.100.kde = sim.kde(10, n=100, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "vm3.100.AICc", ".pdf"))
vm3.100.AICc = sim.AICc(10, n=100, mu=mu3, beta=beta3, pr=pr3, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm3.100.cvISE", ".pdf"))
vm3.100.cvISE = sim.cv(10, n=100, mu=mu3, beta=beta3, pr=pr3, cvinc=5, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm3.100.cvKL", ".pdf"))
vm3.100.cvKL = sim.cv(10, n=100, mu=mu3, beta=beta3, pr=pr3, cvinc=5, 
                      cvloss="KL", seqlen=40, upper=1.09)
dev.off()

# n = 400
pdf(file = paste0(savepath, "vm3.400.kde", ".pdf"))
vm3.400.kde = sim.kde(10, n=400, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "vm3.400.AICc", ".pdf"))
vm3.400.AICc = sim.AICc(10, n=400, mu=mu3, beta=beta3, pr=pr3, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm3.400.cvISE", ".pdf"))
vm3.400.cvISE = sim.cv(10, n=400, mu=mu3, beta=beta3, pr=pr3, cvinc=5, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm3.400.cvKL", ".pdf"))
vm3.400.cvKL = sim.cv(10, n=400, mu=mu3, beta=beta3, pr=pr3, cvinc=5, 
                      cvloss="KL", seqlen=40, upper=1.09)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "vm3.1600.kde", ".pdf"))
vm3.1600.kde = sim.kde(10, n=1600, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "vm3.1600.AICc", ".pdf"))
vm3.1600.AICc = sim.AICc(10, n=1600, mu=mu3, beta=beta3, pr=pr3, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm3.1600.cvISE", ".pdf"))
vm3.1600.cvISE = sim.cv(10, n=1600, mu=mu3, beta=beta3, pr=pr3, cvinc=5, seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "vm3.1600.cvKL", ".pdf"))
vm3.1600.cvKL = sim.cv(10, n=1600, mu=mu3, beta=beta3, pr=pr3, cvinc=5,
                       cvloss="KL", seqlen=40, upper=1.09)
dev.off()

# ----



# wrapped normal mixture --------------------------------------------------

### 1 component
n1 = 400; mu1 = pi * 0; beta1 = 0.8; pr1 = 1
chist(rnpwn(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwn(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

x = rnpwn(100, mu1, beta1, pr1)
fbw(x, 40, 1.07)

# n = 100
pdf(file = paste0(savepath, "wn1.100.kde", ".pdf"))
wn1.100.kde = sim.kde(10, simfam = "npwn", n=100, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "wn1.100.AICc", ".pdf"))
wn1.100.AICc = sim.AICc(10, simfam = "npwn", n=100, mu=mu1, beta=beta1, 
                        seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn1.100.cvISE", ".pdf"))
wn1.100.cvISE = sim.cv(10, simfam = "npwn", n=100, mu=mu1, beta=beta1, 
                       seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn1.100.cvKL", ".pdf"))
wn1.100.cvKL = sim.cv(10, simfam = "npwn", n=100, mu=mu1, beta=beta1, 
                      cvloss="KL", seqlen=40, upper=1.07)
dev.off()

# n = 400
pdf(file = paste0(savepath, "wn1.400.kde", ".pdf"))
wn1.400.kde = sim.kde(10, simfam = "npwn", n=400, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "wn1.400.AICc", ".pdf"))
wn1.400.AICc = sim.AICc(10, simfam = "npwn", n=400, mu=mu1, beta=beta1, 
                        seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn1.400.cvISE", ".pdf"))
wn1.400.cvISE = sim.cv(10, simfam = "npwn", n=400, mu=mu1, beta=beta1, 
                       seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn1.400.cvKL", ".pdf"))
wn1.400.cvKL = sim.cv(10, simfam = "npwn", n=400, mu=mu1, beta=beta1, 
                      cvloss="KL", seqlen=40, upper=1.07)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "wn1.1600.kde", ".pdf"))
wn1.1600.kde = sim.kde(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "wn1.1600.AICc", ".pdf"))
wn1.1600.AICc = sim.AICc(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1,
                         seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn1.1600.cvISE", ".pdf"))
wn1.1600.cvISE = sim.cv(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1,
                        seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn1.1600.cvKL", ".pdf"))
wn1.1600.cvKL = sim.cv(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1, 
                       cvloss="KL", seqlen=40, upper=1.07)
dev.off()



### 2 component
n2 = 400; mu2 = pi * c(0.2, 0.7); beta2 = 0.9; pr2 = c(1/3, 2/3)
chist(rnpwn(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwn(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

x = rnpwn(100, mu2, beta2, pr2)
fbw(x, 40, 1.09)

# n = 100
pdf(file = paste0(savepath, "wn2.100.kde", ".pdf"))
wn2.100.kde = sim.kde(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "wn2.100.AICc", ".pdf"))
wn2.100.AICc = sim.AICc(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2,
                        seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "wn2.100.cvISE", ".pdf"))
wn2.100.cvISE = sim.cv(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2,
                       seqlen=40, upper=1.09)
dev.off()

pdf(file = paste0(savepath, "wn2.100.cvKL", ".pdf"))
wn2.100.cvKL = sim.cv(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2, 
                      cvloss="KL", seqlen=40, upper=1.09)
dev.off()

# n = 400
pdf(file = paste0(savepath, "wn2.400.kde", ".pdf"))
wn2.400.kde = sim.kde(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "wn2.400.AICc", ".pdf"))
wn2.400.AICc = sim.AICc(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2,
                        seqlen=40, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wn2.400.cvISE", ".pdf"))
wn2.400.cvISE = sim.cv(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2,
                       seqlen=40, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wn2.400.cvKL", ".pdf"))
wn2.400.cvKL = sim.cv(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2, 
                      cvloss="KL", seqlen=40, upper=1.08)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "wn2.1600.kde", ".pdf"))
wn2.1600.kde = sim.kde(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "wn2.1600.AICc", ".pdf"))
wn2.1600.AICc = sim.AICc(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2,
                         seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn2.1600.cvISE", ".pdf"))
wn2.1600.cvISE = sim.cv(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2,
                        seqlen=40, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wn2.1600.cvKL", ".pdf"))
wn2.1600.cvKL = sim.cv(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2, 
                       cvloss="KL", seqlen=40, upper=1.07)
dev.off()



### 3 components
n3 = 400; mu3 = pi * c(0.2, 0.7, 1.1); beta3 = 0.9; pr3 = c(3,1,2)/6
chist(rnpwn(n3, mu3, beta3, pr3), nlabels = 0, nbins = 100)
f3 = function(x) dmixwn(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE, col = 4)

x = rnpwn(100, mu3, beta3, pr3)
fbw(x, 40, 1.11)

# n = 100
pdf(file = paste0(savepath, "wn3.100.kde", ".pdf"))
wn3.100.kde = sim.kde(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "wn3.100.AICc", ".pdf"))
wn3.100.AICc = sim.AICc(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3, 
                        seqlen=40, upper=1.11)
dev.off()

pdf(file = paste0(savepath, "wn3.100.cvISE", ".pdf"))
wn3.100.cvISE = sim.cv(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3, 
                       seqlen=40, upper=1.11)
dev.off()

pdf(file = paste0(savepath, "wn3.100.cvKL", ".pdf"))
wn3.100.cvKL = sim.cv(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3, 
                      cvloss="KL", seqlen=40, upper=1.11)
dev.off()

# n = 400
pdf(file = paste0(savepath, "wn3.400.kde", ".pdf"))
wn3.400.kde = sim.kde(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "wn3.400.AICc", ".pdf"))
wn3.400.AICc = sim.AICc(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3, 
                        seqlen=40, upper=1.11)
dev.off()

pdf(file = paste0(savepath, "wn3.400.cvISE", ".pdf"))
wn3.400.cvISE = sim.cv(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3,
                       seqlen=40, upper=1.11)
dev.off()

pdf(file = paste0(savepath, "wn3.400.cvKL", ".pdf"))
wn3.400.cvKL = sim.cv(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3, 
                      cvloss="KL", seqlen=40, upper=1.11)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "wn3.1600.kde", ".pdf"))
wn3.1600.kde = sim.kde(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "wn3.1600.AICc", ".pdf"))
wn3.1600.AICc = sim.AICc(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3,
                         seqlen=40, upper=1.1)
dev.off()

pdf(file = paste0(savepath, "wn3.1600.cvISE", ".pdf"))
wn3.1600.cvISE = sim.cv(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3,
                        seqlen=40, upper=1.1)
dev.off()

pdf(file = paste0(savepath, "wn3.1600.cvKL", ".pdf"))
wn3.1600.cvKL = sim.cv(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3, 
                       cvloss="KL", seqlen=40, upper=1.1)
dev.off()

# ----



# wrapped cauchy mixture --------------------------------------------------

### 1 component
n1 = 400; mu1 = pi * 0; beta1 = 0.8; pr1 = 1
chist(rnpwc(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwc(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

x = rnpwc(100, mu1, beta1, pr1)
fbw(x, 40, 1.07)

# n = 100
pdf(file = paste0(savepath, "wc1.100.kde", ".pdf"))
wc1.100.kde = sim.kde(10, simfam = "npwc", n=100, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "wc1.100.AICc", ".pdf"))
wc1.100.AICc = sim.AICc(10, simfam = "npwc", n=100, mu=mu1, beta=beta1, 
                        seqlen=40, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wc1.100.cvISE", ".pdf"))
wc1.100.cvISE = sim.cv(10, simfam = "npwc", n=100, mu=mu1, beta=beta1,
                       seqlen=40, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wc1.100.cvKL", ".pdf"))
wc1.100.cvKL = sim.cv(10, simfam = "npwc", n=100, mu=mu1, beta=beta1, 
                      cvloss="KL", seqlen=40, upper=1.06)
dev.off()
Sys.time()

# n = 400
pdf(file = paste0(savepath, "wc1.400.kde", ".pdf"))
wc1.400.kde = sim.kde(10, simfam = "npwc", n=400, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "wc1.400.AICc", ".pdf"))
wc1.400.AICc = sim.AICc(10, simfam = "npwc", n=400, mu=mu1, beta=beta1,
                        seqlen=45, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wc1.400.cvISE", ".pdf"))
wc1.400.cvISE = sim.cv(10, simfam = "npwc", n=400, mu=mu1, beta=beta1,
                       seqlen=45, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wc1.400.cvKL", ".pdf"))
wc1.400.cvKL = sim.cv(10, simfam = "npwc", n=400, mu=mu1, beta=beta1, 
                      cvloss="KL", seqlen=45, upper=1.07)
dev.off()
Sys.time()

# n = 1600
pdf(file = paste0(savepath, "wc1.1600.kde", ".pdf"))
wc1.1600.kde = sim.kde(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1)
dev.off()

pdf(file = paste0(savepath, "wc1.1600.AICc", ".pdf"))
wc1.1600.AICc = sim.AICc(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1, 
                         seqlen=45, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc1.1600.cvISE", ".pdf"))
wc1.1600.cvISE = sim.cv(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1, 
                        seqlen=45, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc1.1600.cvKL", ".pdf"))
wc1.1600.cvKL = sim.cv(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1, 
                       cvloss="KL", seqlen=45, upper=1.08)
dev.off()
Sys.time()



### 2 components
n2 = 400; mu2 = pi * c(0.2, 0.5); beta2 = 0.8; pr2 = 1:2/3
chist(rnpwc(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwc(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

x = rnpwc(100, mu2, beta2, pr2)
fbw(x, 45, 1.06)

# n = 100
pdf(file = paste0(savepath, "wc2.100.kde", ".pdf"))
wc2.100.kde = sim.kde(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "wc2.100.AICc", ".pdf"))
wc2.100.AICc = sim.AICc(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2,
                        seqlen=45, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wc2.100.cvISE", ".pdf"))
wc2.100.cvISE = sim.cv(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2,
                       seqlen=45, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wc2.100.cvKL", ".pdf"))
wc2.100.cvKL = sim.cv(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2, 
                      cvloss="KL", seqlen=45, upper=1.06)
dev.off()
Sys.time()

# n = 400
pdf(file = paste0(savepath, "wc2.400.kde", ".pdf"))
wc2.400.kde = sim.kde(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "wc2.400.AICc", ".pdf"))
wc2.400.AICc = sim.AICc(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2,
                        seqlen=45, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wc2.400.cvISE", ".pdf"))
wc2.400.cvISE = sim.cv(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2,
                       seqlen=45, upper=1.07)
dev.off()

pdf(file = paste0(savepath, "wc2.400.cvKL", ".pdf"))
wc2.400.cvKL = sim.cv(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2, 
                      cvloss="KL", seqlen=45, upper=1.07)
dev.off()
Sys.time()

# n = 1600
pdf(file = paste0(savepath, "wc2.1600.kde", ".pdf"))
wc2.1600.kde = sim.kde(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2)
dev.off()

pdf(file = paste0(savepath, "wc2.1600.AICc", ".pdf"))
wc2.1600.AICc = sim.AICc(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2, 
                         seqlen=45, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc2.1600.cvISE", ".pdf"))
wc2.1600.cvISE = sim.cv(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2, 
                        seqlen=45, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc2.1600.cvKL", ".pdf"))
wc2.1600.cvKL = sim.cv(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2, 
                       cvloss="KL", seqlen=45, upper=1.08)
dev.off()
Sys.time()



### 3 components
n3 = 400; mu3 = pi * c(0.2, 0.7, 1); beta3 = 0.8; pr3 = c(3,1,2)/6
chist(rnpwc(n3, mu3, beta3, pr3), nlabels = 0, nbins = 100)
f3 = function(x) dmixwc(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE, col = 4)

x = rnpwc(100, mu3, beta3, pr3)
fbw(x, 50, 1.08)

# n = 100
pdf(file = paste0(savepath, "wc3.100.kde", ".pdf"))
wc3.100.kde = sim.kde(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "wc3.100.AICc", ".pdf"))
wc3.100.AICc = sim.AICc(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3,
                        seqlen=50, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc3.100.cvISE", ".pdf"))
wc3.100.cvISE = sim.cv(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3,
                       seqlen=50, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc3.100.cvKL", ".pdf"))
wc3.100.cvKL = sim.cv(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3, 
                      cvloss="KL", seqlen=50, upper=1.08)
dev.off()
Sys.time()

# n = 400
pdf(file = paste0(savepath, "wc3.400.kde", ".pdf"))
wc3.400.kde = sim.kde(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "wc3.400.AICc", ".pdf"))
wc3.400.AICc = sim.AICc(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3,
                        seqlen=50, upper=1.085)
dev.off()

pdf(file = paste0(savepath, "wc3.400.cvISE", ".pdf"))
wc3.400.cvISE = sim.cv(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3,
                       seqlen=50, upper=1.085)
dev.off()

pdf(file = paste0(savepath, "wc3.400.cvKL", ".pdf"))
wc3.400.cvKL = sim.cv(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3, 
                      cvloss="KL", seqlen=50, upper=1.085)
dev.off()
Sys.time()

# n = 1600
pdf(file = paste0(savepath, "wc3.1600.kde", ".pdf"))
wc3.1600.kde = sim.kde(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3)
dev.off()

pdf(file = paste0(savepath, "wc3.1600.AICc", ".pdf"))
wc3.1600.AICc = sim.AICc(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3,
                         seqlen=60, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc3.1600.cvISE", ".pdf"))
wc3.1600.cvISE = sim.cv(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3,
                        seqlen=60, upper=1.08)
dev.off()

pdf(file = paste0(savepath, "wc3.1600.cvKL", ".pdf"))
wc3.1600.cvKL = sim.cv(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3, 
                       cvloss="KL", seqlen=60, upper=1.08)
dev.off()
Sys.time()

# ----



# wrapped skewed normal ---------------------------------------------------

### 1 component
n1 = 200; xi1 = pi * 0.5; omega1 = 1; alpha1 = -5
chist(rnpwsn(n1, xi1, omega1, alpha1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwsn(x, xi1, omega1, alpha1, pr1)
cdensity(f1, add = TRUE, col = 4)

x = rnpwsn(100, xi1, omega1, alpha1)
fbw(x, 50, 1.05)

# n = 100
pdf(file = paste0(savepath, "wsn1.100.kde", ".pdf"))
wsn1.100.kde = sim.kde(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
dev.off()

pdf(file = paste0(savepath, "wsn1.100.AICc", ".pdf"))
wsn1.100.AICc = sim.AICc(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, 
                         alpha=alpha1, seqlen=50, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wsn1.100.cvISE", ".pdf"))
wsn1.100.cvISE = sim.cv(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, 
                        alpha=alpha1, seqlen=50, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wsn1.100.cvKL", ".pdf"))
wsn1.100.cvKL = sim.cv(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, 
                       alpha=alpha1, cvloss="KL", seqlen=50, upper=1.06)
dev.off()

# n = 400
pdf(file = paste0(savepath, "wsn1.400.kde", ".pdf"))
wsn1.400.kde = sim.kde(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
dev.off()

pdf(file = paste0(savepath, "wsn1.400.AICc", ".pdf"))
wsn1.400.AICc = sim.AICc(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, 
                         alpha=alpha1, seqlen=50, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wsn1.400.cvISE", ".pdf"))
wsn1.400.cvISE = sim.cv(10, n=400, simfam="npwsn", mu=xi1,beta=omega1, 
                        alpha=alpha1, seqlen=50, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wsn1.400.cvKL", ".pdf"))
wsn1.400.cvKL = sim.cv(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, 
                       alpha=alpha1, cvloss="KL", seqlen=50, upper=1.06)
dev.off()

# n = 1600
pdf(file = paste0(savepath, "wsn1.1600.kde", ".pdf"))
wsn1.1600.kde = sim.kde(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
dev.off()

pdf(file = paste0(savepath, "wsn1.1600.AICc", ".pdf"))
wsn1.1600.AICc = sim.AICc(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, 
                          alpha=alpha1, seqlen=50, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wsn1.1600.cvISE", ".pdf"))
wsn1.1600.cvISE = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, 
                         alpha=alpha1, seqlen=50, upper=1.06)
dev.off()

pdf(file = paste0(savepath, "wsn1.1600.cvKL", ".pdf"))
wsn1.1600.cvKL = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, 
                        alpha=alpha1, cvloss="KL", seqlen=50, upper=1.06)
dev.off()

# ----



# real data ---------------------------------------------------------------

?Turtles_radians
?HurricanesGulfofMexico1971to2008
?ProteinsAAA

# CircNNTSR::turtle
data(Turtles_radians)
turtle = Turtles_radians
length(turtle)
chist(turtle, nbins = 100, nlabels = 0)

fbw(turtle, seqlen = 40, upper = 1.05)

pdf(file = paste0(savepath2, "turtle.kde", ".pdf"))
turtle.kde = cv.real.kde(turtle, repetition=10)
dev.off()

pdf(file = paste0(savepath2, "turtle.AICc", ".pdf"))
turtle.AICc = cv.real.mde(turtle, "AICc", repetition=10, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath2, "turtle.cv", ".pdf"))
turtle.cv = cv.real.mde(turtle, "cv", repetition=10, cvinc=0.1, seqlen=40, upper=1.05)
dev.off()

pdf(file = paste0(savepath2, "turtle.cv.KL", ".pdf"))
turtle.cv.KL = cv.real.mde(turtle, "cv", repetition=10, cvinc=0.1, 
                           cvloss="KL", seqlen=40, upper=1.05)
dev.off()
Sys.time()



# CircNNTSR hurricane2
data(HurricanesGulfofMexico1971to2008)
hurricane = HurricanesGulfofMexico1971to2008$V1 * 2 * pi
length(hurricane)
chist(hurricane, nbins = 100)

fbw(hurricane, seqlen = 50, upper = 1.05)

pdf(file = paste0(savepath2, "hurricane.kde", ".pdf"))
hurricane.kde = cv.real.kde(hurricane, repetition=10)
dev.off()

pdf(file = paste0(savepath2, "hurricane.AICc", ".pdf"))
hurricane.AICc = cv.real.mde(hurricane, "AICc", repetition=10, seqlen=50, upper=1.05)
dev.off()

pdf(file = paste0(savepath2, "hurricane.cv", ".pdf"))
hurricane.cv = cv.real.mde(hurricane, "cv", repetition=10, seqlen=50, upper=1.05)
dev.off()

pdf(file = paste0(savepath2, "hurricane.cv.KL", ".pdf"))
hurricane.cv.KL = cv.real.mde(hurricane, "cv", repetition=10, cvloss="KL", 
                              seqlen=50, upper=1.05)
dev.off()



# CircNNTSR protein
data(ProteinsAAA)
dim(ProteinsAAA)
# protein1 = ProteinsAAA$V1
# chist(protein1,  nbins = 100, nlabels = 0)
protein = ProteinsAAA$V2
chist(protein, nbins = 100, nlabels = 0)

fbw(protein, seqlen = 100, upper = 1.045)

pdf(file = paste0(savepath2, "protein.kde", ".pdf"))
protein.kde = cv.real.kde(protein, repetition=10)
dev.off()

pdf(file = paste0(savepath2, "protein.AICc", ".pdf"))
protein.AICc = cv.real.mde(protein, "AICc", repetition=10, seqlen=100, upper=1.045)
dev.off()

pdf(file = paste0(savepath2, "protein.cv", ".pdf"))
protein.cv = cv.real.mde(protein, "cv", repetition=10, cvinc=0.1, 
                         seqlen=100, upper=1.045)
dev.off()

pdf(file = paste0(savepath2, "protein.cv.KL", ".pdf"))
protein.cv.KL = cv.real.mde(protein, "cv", repetition=10, cvinc=0.1, 
                            cvloss="KL", seqlen=100, upper=1.045)
dev.off()



# dragonfly
data(dragonfly)
dfly = dragonfly$orientation
chist(dfly, nbins = 100, nlabels = 0)

fbw(dfly, seqlen = 80, upper = 1.08)

pdf(file = paste0(savepath2, "dfly.kde", ".pdf"))
dfly.kde = cv.real.kde(dfly, repetition=10)
dev.off()

pdf(file = paste0(savepath2, "dfly.AICc", ".pdf"))
dfly.AICc = cv.real.mde(dfly, "AICc", repetition=10, seqlen=80, upper=1.08)
dev.off()

pdf(file = paste0(savepath2, "dfly.cv", ".pdf"))
dfly.cv = cv.real.mde(dfly, "cv", repetition=10, cvinc=0.1, seqlen=80, upper=1.08)
dev.off()

pdf(file = paste0(savepath2, "dfly.cv.KL", ".pdf"))
dfly.cv.KL = cv.real.mde(dfly, "cv", repetition=10, cvinc=0.1, cvloss="KL", 
                         seqlen=80, upper=1.08)
dev.off()



# cross bed
crossbed = unlist(fisherB6) / 180 * pi
length(crossbed)
chist(crossbed, nbins = 100, nlabels = 0)

fbw(crossbed, seqlen = 45, upper = 1.07)

pdf(file = paste0(savepath2, "crossbed.kde", ".pdf"))
crossbed.kde = cv.real.kde(crossbed, repetition=10)
dev.off()

pdf(file = paste0(savepath2, "crossbed.AICc", ".pdf"))
crossbed.AICc = cv.real.mde(crossbed, "AICc", repetition=10, seqlen=45, upper=1.07)
dev.off()

pdf(file = paste0(savepath2, "crossbed.cv", ".pdf"))
crossbed.cv = cv.real.mde(crossbed, "cv", repetition=10, cvinc=0.1, 
                          seqlen=45, upper=1.07)
dev.off()

pdf(file = paste0(savepath2, "crossbed.cv.KL", ".pdf"))
crossbed.cv.KL = cv.real.mde(crossbed, "cv", repetition=10, cvinc=0.1, 
                             cvloss="KL", seqlen=45, upper=1.07)
dev.off()
Sys.time()



### icu
data(fisherB1c)
icu = unclass(fisherB1c)
attributes(icu) = NULL
icu = icu / 24 * 2 * pi

fbw(icu, seqlen = 60, upper = 1.09)

pdf(file = paste0(savepath2, "icu.kde", ".pdf"))
icu.kde = cv.real.kde(icu, repetition=10)
dev.off()

pdf(file = paste0(savepath2, "icu.AICc", ".pdf"))
icu.AICc = cv.real.mde(icu, "AICc", repetition=10, seqlen=60, upper=1.09)
dev.off()

pdf(file = paste0(savepath2, "icu.cv", ".pdf"))
icu.cv = cv.real.mde(icu, "cv", repetition=10, cvinc=0.1, seqlen=60, upper=1.095)
dev.off()

pdf(file = paste0(savepath2, "icu.cv.KL", ".pdf"))
icu.cv.KL = cv.real.mde(icu, "cv", repetition=10, cvinc=0.1, cvloss="KL", 
                        seqlen=60, upper=1.09)
dev.off()
Sys.time()



# ### ant
# data(Ants_radians)
# ant = Ants_radians
# 
# pdf(file = paste0(savepath2, "ant.kde", ".pdf"))
# ant.kde = cv.real.kde(ant, repetition=10)
# dev.off()
# 
# pdf(file = paste0(savepath2, "ant.AICc", ".pdf"))
# ant.AICc = cv.real.mde(ant, "AICc", repetition=10, seqlen = 20, upper = 1.5)
# dev.off()
# 
# pdf(file = paste0(savepath2, "ant.cv", ".pdf"))
# ant.cv = cv.real.mde(ant, "cv", repetition=10, cvinc = 0.1, seqlen = 20, upper = 1.5)
# dev.off()
# 
# pdf(file = paste0(savepath2, "ant.cv.KL", ".pdf"))
# ant.cv.KL = cv.real.mde(ant, "cv", repetition=10, cvloss = "KL", cvinc = 0.1, seqlen = 20, upper = 1.5)
# dev.off()
# Sys.time()



# # wind
# data(wind, package = "NPCirc")
# length(wind$wind.dir)
# winddir = wind$wind.dir
# 
# pdf(file = paste0(savepath2, "winddir.kde", ".pdf"))
# winddir.kde = cv.real.kde(winddir, repetition=10)
# dev.off()
# 
# pdf(file = paste0(savepath2, "winddir.AICc", ".pdf"))
# winddir.AICc = cv.real.mde(winddir, "AICc", repetition=10, seqlen = 30, upper = 1.25)
# dev.off()
# 
# pdf(file = paste0(savepath2, "winddir.cv", ".pdf"))
# winddir.cv = cv.real.mde(winddir, "cv", repetition=10, cvinc = 0.1, 
#                          seqlen = 26, upper = 1.25)
# dev.off()
# 
# pdf(file = paste0(savepath2, "winddir.cv.KL", ".pdf"))
# winddir.cv.KL = cv.real.mde(winddir, "cv", repetition=10, cvloss = "KL", 
#                             cvinc = 0.1, seqlen = 26, upper = 1.25)
# dev.off()
# Sys.time()

# ----



# debug -------------------------------------------------------------------

saveRDS({rnorm(10); "a" = print(1)}, file = paste0(savepath3, "rnorm.rds"))
readRDS(file = paste0(savepath3, "rnorm.rds"))

filename = paste0(savepath3, "wsn1.1600.cvISE", ".rds")
saveRDS(wsn1.1600.cvISE = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1), 
        file = filename)
# readRDS(file = filename)

filename = paste0(savepath3, "wsn1.1600.cvKL", ".rds")
saveRDS(wsn1.1600.cvKL = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1, cvloss="KL"), 
        file = filename)

# readRDS(file = filename)



pdf(file = paste0(savepath, "wsn1.1600.cvISE", ".pdf"))
wsn1.1600.cvISE = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
dev.off()

pdf(file = paste0(savepath, "wsn1.1600.cvKL", ".pdf"))
wsn1.1600.cvKL = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1, cvloss="KL")
dev.off()

filename = paste0(savepath3, "wsn1.cvISE.10", ".rds")
saveRDS({pdf(file = paste0(savepath3, "wsn1.cvISE.10", ".pdf")); 
  wsn1.cvISE.10 = sim.cv(10, n=10, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1);
  dev.off()}, 
  file = filename)
readRDS(file = filename)

# pdf then saveRDS
pdf(file = paste0(savepath3, "wsn1.cvISE.10", ".pdf"))
filename = paste0(savepath3, "wsn1.cvISE.10", ".rds")
saveRDS(sim.cv(10, n=10, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1), 
        file = filename)
dev.off()
readRDS(file = filename)

# ----



# results -----------------------------------------------------------------

getres("vm1", 100)
getres("vm1", 400)
getres("vm1", 1600)

getres("vm2", 100)
getres("vm2", 400)
getres("vm2", 1600)

getres("vm3", 100)
getres("vm3", 400)
getres("vm3", 1600)

getres("wn1", 100)
getres("wn1", 400)
getres("wn1", 1600)

getres("wn2", 100)
getres("wn2", 400)
getres("wn2", 1600)

getres("wn3", 100)
getres("wn3", 400)
getres("wn3", 1600)

getres("wc1", 100)
getres("wc1", 400)
getres("wc1", 1600)

getres("wc2", 100)
getres("wc2", 400)
getres("wc2", 1600)

getres("wc3", 100)
getres("wc3", 400)
getres("wc3", 1600)

getres("wsn1", 100)
getres("wsn1", 400)
getres("wsn1", 1600)

getres("wsn2", 100)
getres("wsn2", 400)
getres("wsn2", 1600)

getres("wsn3", 100)
getres("wsn3", 400)
getres("wsn3", 1600)



# aggregate results
rbind(getres("vm1", 100), getres("vm1", 400), getres("vm1", 1600))

array(c(getres("vm1", 100), getres("vm1", 400), getres("vm1", 1600)),
      dim = c(7, 3, 3), 
      dimnames = list("method" = 1:7, "loss" = c("ISE", "KL", "HD"), 
                      "size" = c(100, 400, 1600)))

# not using
aggres = function(family="vm1", loss=c("ISE","KL","HD"),  
                  type=c("mean","se","numcomp","betavalue"), 
                  mul=100, digit=2) {
  type = match.arg(type)
  loss = match.arg(loss)
  numloss = switch(loss, "ISE" = 1, "KL" = 2, "HD" = 3)
  
  m100 = getres(family=family, size=100, type=type, mul=mul, digit=digit)
  m400 = getres(family=family, size=400, type=type, mul=mul, digit=digit)
  m1600 = getres(family=family, size=1600, type=type, mul=mul, digit=digit)
  m100 = m100[c(7:4, 1:3), ]
  m400 = m400[c(7:4, 1:3), ]
  m1600 = m1600[c(7:4, 1:3), ]
  
  ar = array(c(m100, m400, m1600),
             dim = c(7, 3, 3), 
             dimnames = list("method" = c("RT","PI","LSCV","LCV",
                                          "AICc","cvISE","cvKL"), 
                             "loss" = c("ISE", "KL", "HD"), 
                             "size" = c(100, 400, 1600)))
  ar[, numloss, ]
}

# aggres
aggres = function(family="vm1", type=c("mean","se","numcomp","betavalue"), 
                  mul=100, digit=2) {
  type = match.arg(type)

  m100 = getres(family=family, size=100, type=type, mul=mul, digit=digit)
  m400 = getres(family=family, size=400, type=type, mul=mul, digit=digit)
  m1600 = getres(family=family, size=1600, type=type, mul=mul, digit=digit)
  m100 = m100[c(7:4, 1:3), ]
  m400 = m400[c(7:4, 1:3), ]
  m1600 = m1600[c(7:4, 1:3), ]
  
  array(c(m100, m400, m1600),
        dim = c(7, 3, 3), 
        dimnames = list("method" = c("RT","PI","LSCV","LCV","AICc","cvISE","cvKL"), 
                        "loss" = c("ISE", "KL", "HD"),
                        "size" = c(100, 400, 1600)))
}

aggres("vm1")
aggres("vm3")
aggres("wn2")
aggres("wc3")
aggres("wsn1")

aggres("wsn1", "se")[, 1, ]

# convert to table in latex
fmin = function(x) which(x == min(x))

latexres = function(family="vm1", loss=c("ISE","KL","HD"), mul=100, digit=2) {
  loss = match.arg(loss)
  nloss = switch(loss, "ISE" = 1, "KL" = 2, "HD" = 3)
  
  mat = NULL
  for (i in 1:length(family)) {
    mean = aggres(family = family[i], type = "mean", mul = mul, digit = digit)[, nloss, ]
    se = aggres(family = family[i], type = "se", mul = mul, digit = digit)[, nloss, ]
    new = paste(format(c(mean), digits = digit),
                paste0("(", format(se, digits = digit), ")"), sep = " ")
    # ind = apply(mean, 2, which.min) + 7 * 0:2  # only one minimum
    ind = unlist(mapply("+", apply(mean, 2, fmin), 7 * 0:2))  # find multiple minimum
    new[ind] = paste0("\\textbf{", new[ind], "}")
    mat = cbind(mat, new)
  }
  c1 = c("KDE", rep("", 3), "MDE", rep("", 2))
  c2 = c("RT", "PI", "LSCV", "LCV", "AICc", "CVISE", "CVKL")
  mat = cbind(c1, c2, mat)
  result = apply(mat, 1, paste, collapse = " & ")
  result = paste0(result, "\\\\")
  res = character(length(result) + 10)
  res[-c(1:2, 7, 11:12, 17, 21:22, 27, 31)] = result
  res[1] = "\\toprule \\textit{Estimator} & \\textit{Selection method} & vM2 & WN3 & WC2 & WSN \\\\ \\toprule"
  res[c(2, 12, 22)] = paste0("&&\\multicolumn{4}{c}{$n = ", 
                             c(100,400,1600), 
                             "$} \\\\ \\cmidrule{3-6}")
  res[c(7, 17, 27)] = "\\midrule"
  res[c(11, 21, 31)] = "\\bottomrule"
  cat(res, fill = TRUE)
}

latexres(c("vm2", "wn3", "wc2", "wsn1"))
latexres(c("vm2", "wn3", "wc2", "wsn1"), loss = "KL")
latexres(c("vm2", "wn3", "wc2", "wsn1"), loss = "HD")

# ----







