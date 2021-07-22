
# von Mises mixture -------------------------------------------------------

### 1 component
n1 = 400; mu1 = pi * 0; beta1 = 20; pr1 = 1
chist(rnpvm(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixvm(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

# n = 100
vm1.100.kde = sim.kde(10, n=100, mu=mu1, beta=beta1)
vm1.100.AICc = sim.AICc(10, n=100, mu=mu1, beta=beta1)
vm1.100.cvISE = sim.cv(10, n=100, mu=mu1, beta=beta1)
vm1.100.cvKL = sim.cv(10, n=100, mu=mu1, beta=beta1, cvloss="KL")

getres("vm1", 100)
getres("vm1", 100, "se")
getres("vm1", 100, "numcomp")

# n = 400
vm1.400.kde = sim.kde(10, n=400, mu=mu1, beta=beta1)
vm1.400.AICc = sim.AICc(10, n=400, mu=mu1, beta=beta1)
vm1.400.cvISE = sim.cv(10, n=400, mu=mu1, beta=beta1)
vm1.400.cvKL = sim.cv(10, n=400, mu=mu1, beta=beta1, cvloss="KL")

getres("vm1", 400)
getres("vm1", 400, "se")
getres("vm1", 400, "numcomp")

# n = 1600
vm1.1600.kde = sim.kde(10, n=1600, mu=mu1, beta=beta1)
vm1.1600.AICc = sim.AICc(10, n=1600, mu=mu1, beta=beta1)
vm1.1600.cvISE = sim.cv(10, n=1600, mu=mu1, beta=beta1)
vm1.1600.cvKL = sim.cv(10, n=1600, mu=mu1, beta=beta1, cvloss="KL")

getres("vm1", 1600)
getres("vm1", 1600, "se")
getres("vm1", 1600, "numcomp")



## 2 components
n2 = 400; mu2 = pi * c(0.2, 0.5); beta2 = 15; pr2 = c(1, 2)/3
chist(rnpvm(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixvm(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

# n = 100
vm2.100.kde = sim.kde(10, n=100, mu=mu2, beta=beta2, pr=pr2)
vm2.100.AICc = sim.AICc(10, n=100, mu=mu2, beta=beta2, pr=pr2)
vm2.100.cvISE = sim.cv(10, n=100, mu=mu2, beta=beta2, pr=pr2)
vm2.100.cvKL = sim.cv(10, n=100, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("vm2", 100)
getres("vm2", 100, "se")
getres("vm2", 100, "numcomp")

# n = 400
vm2.400.kde = sim.kde(10, n=400, mu=mu2, beta=beta2, pr=pr2)
vm2.400.AICc = sim.AICc(10, n=400, mu=mu2, beta=beta2, pr=pr2)
vm2.400.cvISE = sim.cv(10, n=400, mu=mu2, beta=beta2, pr=pr2)
vm2.400.cvKL = sim.cv(10, n=400, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("vm2", 400)
getres("vm2", 400, "se")
getres("vm2", 400, "numcomp")

# n = 1600
vm2.1600.kde = sim.kde(10, n=1600, mu=mu2, beta=beta2, pr=pr2)
vm2.1600.AICc = sim.AICc(10, n=1600, mu=mu2, beta=beta2, pr=pr2)
vm2.1600.cvISE = sim.cv(10, n=1600, mu=mu2, beta=beta2, pr=pr2)
vm2.1600.cvKL = sim.cv(10, n=1600, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("vm2", 1600)
getres("vm2", 1600, "se")
getres("vm2", 1600, "numcomp")



### 3 components
n3 = 400; mu3 = pi * c(0.2, 0.5, 0.8); beta3 = 15; pr3 = c(1,3,2)/6
chist(rnpvm(n3, mu3, beta3, pr3), nlabels = 0, nbins = 100)
f3 = function(x) dmixvm(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE, col = 4)

# n = 100
vm3.100.kde = sim.kde(10, n=100, mu=mu3, beta=beta3, pr=pr3)
vm3.100.AICc = sim.AICc(10, n=100, mu=mu3, beta=beta3, pr=pr3)
vm3.100.cvISE = sim.cv(10, n=100, mu=mu3, beta=beta3, pr=pr3, cvinc=5)
vm3.100.cvKL = sim.cv(10, n=100, mu=mu3, beta=beta3, pr=pr3, cvinc=5, cvloss="KL")

getres("vm3", 100)
getres("vm3", 100, "se")
getres("vm3", 100, "numcomp")

# n = 400
vm3.400.kde = sim.kde(10, n=400, mu=mu3, beta=beta3, pr=pr3)
vm3.400.AICc = sim.AICc(10, n=400, mu=mu3, beta=beta3, pr=pr3)
vm3.400.cvISE = sim.cv(10, n=400, mu=mu3, beta=beta3, pr=pr3, cvinc=5)
vm3.400.cvKL = sim.cv(10, n=400, mu=mu3, beta=beta3, pr=pr3, cvinc=5, cvloss="KL")

getres("vm3", 400)
getres("vm3", 400, "se")
getres("vm3", 400, "numcomp")

# n = 1600
vm3.1600.kde = sim.kde(10, n=1600, mu=mu3, beta=beta3, pr=pr3)
vm3.1600.AICc = sim.AICc(10, n=1600, mu=mu3, beta=beta3, pr=pr3)
vm3.1600.cvISE = sim.cv(10, n=1600, mu=mu3, beta=beta3, pr=pr3, cvinc=5)
vm3.1600.cvKL = sim.cv(10, n=1600, mu=mu3, beta=beta3, pr=pr3, cvinc=5, cvloss="KL")

getres("vm3", 1600)
getres("vm3", 1600, "se")
getres("vm3", 1600, "numcomp")

# ----



# wrapped normal mixture --------------------------------------------------

### 1 component
n1 = 400; mu1 = pi * 0; beta1 = 0.8; pr1 = 1
chist(rnpwn(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwn(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

# n = 100
wn1.100.kde = sim.kde(10, simfam = "npwn", n=100, mu=mu1, beta=beta1)
wn1.100.AICc = sim.AICc(10, simfam = "npwn", n=100, mu=mu1, beta=beta1)
wn1.100.cvISE = sim.cv(10, simfam = "npwn", n=100, mu=mu1, beta=beta1)
wn1.100.cvKL = sim.cv(10, simfam = "npwn", n=100, mu=mu1, beta=beta1, cvloss="KL")

getres("wn1", 100)
getres("wn1", 100, "se")
getres("wn1", 100, "numcomp")
getres("wn1", 100, "betavalue")

# n = 400
wn1.400.kde = sim.kde(10, simfam = "npwn", n=400, mu=mu1, beta=beta1)
wn1.400.AICc = sim.AICc(10, simfam = "npwn", n=400, mu=mu1, beta=beta1)
wn1.400.cvISE = sim.cv(10, simfam = "npwn", n=400, mu=mu1, beta=beta1)
wn1.400.cvKL = sim.cv(10, simfam = "npwn", n=400, mu=mu1, beta=beta1, cvloss="KL")

getres("wn1", 400)
getres("wn1", 400, "se")
getres("wn1", 400, "betavalue")

# n = 1600
wn1.1600.kde = sim.kde(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1)
wn1.1600.AICc = sim.AICc(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1)
wn1.1600.cvISE = sim.cv(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1)
wn1.1600.cvKL = sim.cv(10, simfam = "npwn", n=1600, mu=mu1, beta=beta1, cvloss="KL")

getres("wn1", 1600, mul = 1000)
getres("wn1", 1600, "se", mul = 1000)
getres("wn1", 1600, "betavalue")



### 2 component
n2 = 400; mu2 = pi * c(0.2, 0.7); beta2 = 0.9; pr2 = c(1/3, 2/3)
# n2 = 400; mu2 = pi * c(0.2, 2); beta2 = 0.9; pr2 = c(1/3, 2/3)
# n2 = 400; mu2 = pi * c(0.2, 0.6); beta2 = 0.9; pr2 = c(1/3, 2/3)
chist(rnpwn(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwn(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

# n = 100
wn2.100.kde = sim.kde(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2)
wn2.100.AICc = sim.AICc(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2)
wn2.100.cvISE = sim.cv(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2)
wn2.100.cvKL = sim.cv(10, simfam="npwn", n=100, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("wn2", 100)
getres("wn2", 100, "se")
getres("wn2", 100, "betavalue")

# n = 400
wn2.400.kde = sim.kde(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2)
wn2.400.AICc = sim.AICc(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2)
wn2.400.cvISE = sim.cv(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2)
wn2.400.cvKL = sim.cv(10, simfam="npwn", n=400, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("wn2", 400)
getres("wn2", 400, "se")
getres("wn2", 400, "betavalue")

# n = 1600
wn2.1600.kde = sim.kde(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2)
wn2.1600.AICc = sim.AICc(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2)
wn2.1600.cvISE = sim.cv(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2)
wn2.1600.cvKL = sim.cv(10, simfam="npwn", n=1600, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("wn2", 1600)
getres("wn2", 1600, "se")
getres("wn2", 1600, "betavalue")



### 3 components
n3 = 400; mu3 = pi * c(0.2, 0.7, 1.1); beta3 = 0.9; pr3 = c(3,1,2)/6
chist(rnpwn(n3, mu3, beta3, pr3), nlabels = 0, nbins = 100)
f3 = function(x) dmixwn(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE, col = 4)

# n = 100
wn3.100.kde = sim.kde(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3)
wn3.100.AICc = sim.AICc(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3)
wn3.100.cvISE = sim.cv(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3)
wn3.100.cvKL = sim.cv(10, simfam="npwn", n=100, mu=mu3, beta=beta3, pr=pr3, cvloss="KL")

getres("wn3", 100)
getres("wn3", 100, "se")
getres("wn3", 100, "betavalue")

# n = 400
wn3.400.kde = sim.kde(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3)
wn3.400.AICc = sim.AICc(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3)
wn3.400.cvISE = sim.cv(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3)
wn3.400.cvKL = sim.cv(10, simfam="npwn", n=400, mu=mu3, beta=beta3, pr=pr3, cvloss="KL")

getres("wn3", 400)
getres("wn3", 400, "se")
getres("wn3", 400, "betavalue")

# n = 1600
wn3.1600.kde = sim.kde(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3)
wn3.1600.AICc = sim.AICc(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3)
wn3.1600.cvISE = sim.cv(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3)
wn3.1600.cvKL = sim.cv(10, simfam="npwn", n=1600, mu=mu3, beta=beta3, pr=pr3, cvloss="KL")

getres("wn3", 1600)
getres("wn3", 1600, "se")
getres("wn3", 1600, "betavalue")

# ----



# wrapped Cauchy mixture --------------------------------------------------

### 1 component
n1 = 100; mu1 = pi * 0; beta1 = 0.8; pr1 = 1
chist(rnpwc(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwc(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

# n = 100
wc1.100.kde = sim.kde(10, simfam = "npwc", n=100, mu=mu1, beta=beta1)
wc1.100.AICc = sim.AICc(10, simfam = "npwc", n=100, mu=mu1, beta=beta1)
wc1.100.cvISE = sim.cv(10, simfam = "npwc", n=100, mu=mu1, beta=beta1)
wc1.100.cvKL = sim.cv(10, simfam = "npwc", n=100, mu=mu1, beta=beta1, cvloss="KL")

getres("wc1", 100)
getres("wc1", 100, "se")
getres("wc1", 100, "numcomp")
getres("wc1", 100, "betavalue")

# n = 400
wc1.400.kde = sim.kde(10, simfam = "npwc", n=400, mu=mu1, beta=beta1)
wc1.400.AICc = sim.AICc(10, simfam = "npwc", n=400, mu=mu1, beta=beta1)
wc1.400.cvISE = sim.cv(10, simfam = "npwc", n=400, mu=mu1, beta=beta1)
wc1.400.cvKL = sim.cv(10, simfam = "npwc", n=400, mu=mu1, beta=beta1, cvloss="KL")

getres("wc1", 400)
getres("wc1", 400, "se")
getres("wc1", 400, "numcomp")
getres("wc1", 400, "betavalue")

# n = 1600
wc1.1600.kde = sim.kde(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1)
wc1.1600.AICc = sim.AICc(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1)
wc1.1600.cvISE = sim.cv(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1)
wc1.1600.cvKL = sim.cv(10, simfam = "npwc", n=1600, mu=mu1, beta=beta1, cvloss="KL")

getres("wc1", 1600)
getres("wc1", 1600, "se")
getres("wc1", 1600, "numcomp")
getres("wc1", 1600, "betavalue")

set.seed(1789)
q = rnpwc(n=1600, mu=mu1, beta=beta1)
fseq(finvA(q), upper = 1.2, seqlen = 21)
res = cnm(fnpvm(q), init = list(beta = 72))
plot.npvm(fnpvm(q), res$mix, res$beta, dimen = 2)


# test larger beta (smaller sd)
n1 = 500; mu1 = pi * 0; beta1 = 0.9; pr1 = 1
chist(rnpwc(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwc(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

k11 = sim.kde(10, simfam = "npwc", n=500, mu=mu1, beta=beta1)
a11 = sim.AICc(10, simfam = "npwc", n=500, mu=mu1, beta=beta1)
round(rbind(AICc = a11$mean, k11$mean) * 100, 2)

# test smaller beta (larger sd)
n1 = 500; mu1 = pi * 0; beta1 = 0.7; pr1 = 1
chist(rnpwc(n1, mu1, beta1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwc(x, mu1, beta1, pr1)
cdensity(f1, add = TRUE, col = 4)

k12 = sim.kde(10, simfam = "npwc", n=500, mu=mu1, beta=beta1)
a12 = sim.AICc(10, simfam = "npwc", n=500, mu=mu1, beta=beta1)
round(rbind(AICc = a12$mean, k12$mean) * 100, 2)



### 2 components
n2 = 500; mu2 = pi * c(0.2, 0.5); beta2 = 0.8; pr2 = 1:2/3
# n2 = 500; mu2 = pi * c(0.2, 0.8); beta2 = 0.7; pr2 = 1:2/3
chist(rnpwc(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwc(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

# n = 100
wc2.100.kde = sim.kde(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2)
wc2.100.AICc = sim.AICc(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2)
wc2.100.cvISE = sim.cv(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2)
wc2.100.cvKL = sim.cv(10, simfam = "npwc", n=100, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("wc2", 100)
getres("wc2", 100, "se")
getres("wc2", 100, "numcomp")
getres("wc2", 100, "betavalue")

# n = 400
wc2.400.kde = sim.kde(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2)
wc2.400.AICc = sim.AICc(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2)
wc2.400.cvISE = sim.cv(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2)
wc2.400.cvKL = sim.cv(10, simfam = "npwc", n=400, mu=mu2, beta=beta2, pr=pr2, cvloss="KL")

getres("wc2", 400)
getres("wc2", 400, "se")
getres("wc2", 400, "numcomp")

# n = 1600
wc2.1600.kde = sim.kde(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2)
wc2.1600.AICc = sim.AICc(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2, seqlen = 23)
wc2.1600.cvISE = sim.cv(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2, seqlen = 23)
wc2.1600.cvKL = sim.cv(10, simfam = "npwc", n=1600, mu=mu2, beta=beta2, pr=pr2, 
                       cvloss='KL', seqlen = 23)

getres("wc2", 1600)
getres("wc2", 1600, "se")
getres("wc2", 1600, "numcomp")

set.seed(1789)
q = rnpwc(n=1600, mu=mu2, beta=beta2, pr = pr2)
fseq(finvA(q), upper = 1.2, seqlen = 23)
res = cnm(fnpvm(q), init = list(beta = 80))
plot.npvm(fnpvm(q), res$mix, res$beta, dimen = 2)

# test larger beta
n2 = 500; mu2 = pi * c(0.2, 0.5); beta2 = 0.9; pr2 = 1:2/3
chist(rnpwc(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwc(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

k21 = sim.kde(10, simfam = "npwc", n=500, mu=mu2, beta=beta2, pr=pr2)
a21 = sim.AICc(10, simfam = "npwc", n=500, mu=mu2, beta=beta2, pr=pr2)
round(rbind(AICc = a21$mean, k21$mean) * 100, 2)

# test smaller beta
n2 = 500; mu2 = pi * c(0.2, 0.6); beta2 = 0.7; pr2 = 1:2/3
chist(rnpwc(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwc(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

k22 = sim.kde(10, simfam = "npwc", n=500, mu=mu2, beta=beta2, pr=pr2)
a22 = sim.AICc(10, simfam = "npwc", n=500, mu=mu2, beta=beta2, pr=pr2)
round(rbind(AICc = a22$mean, k22$mean) * 100, 2)

# test larger beta
n2 = 500; mu2 = pi * c(0.2, 1.2); beta2 = 0.9; pr2 = 1:2/3
chist(rnpwc(n2, mu2, beta2, pr2), nlabels = 0, nbins = 100)
f2 = function(x) dmixwc(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE, col = 4)

k23 = sim.kde(10, simfam = "npwc", n=500, mu=mu2, beta=beta2, pr=pr2)
a23 = sim.AICc(10, simfam = "npwc", n=500, mu=mu2, beta=beta2, pr=pr2)
round(rbind(AICc = a23$mean, k23$mean) * 100, 2)



### 3 components
n3 = 500; mu3 = pi * c(0.2, 0.7, 1); beta3 = 0.8; pr3 = c(3,1,2)/6
chist(rnpwc(n3, mu3, beta3, pr3), nlabels = 0, nbins = 100)
f3 = function(x) dmixwc(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE, col = 4)

# n = 100
wc3.100.kde = sim.kde(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3)
wc3.100.AICc = sim.AICc(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3)
wc3.100.cvISE = sim.cv(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3)
wc3.100.cvKL = sim.cv(10, simfam = "npwc", n=100, mu=mu3, beta=beta3, pr=pr3, cvloss="KL")

getres("wc3", 100)
getres("wc3", 100, "se")
getres("wc3", 100, "numcomp")

# n = 400
wc3.400.kde = sim.kde(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3)
wc3.400.AICc = sim.AICc(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3)
wc3.400.cvISE = sim.cv(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3)
wc3.400.cvKL = sim.cv(10, simfam = "npwc", n=400, mu=mu3, beta=beta3, pr=pr3, cvloss="KL")

getres("wc3", 400)
getres("wc3", 400, "se")
getres("wc3", 400, "numcomp")

# n = 1600
wc3.1600.kde = sim.kde(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3)
wc3.1600.AICc = sim.AICc(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3)
wc31600.cvISE = sim.cv(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3)
wc3.1600.cvKL = sim.cv(10, simfam = "npwc", n=1600, mu=mu3, beta=beta3, pr=pr3, cvloss="KL")

getres("wc3", 1600)
getres("wc3", 1600, "se")
getres("wc3", 1600, "numcomp")

set.seed(1789)
q = rnpwc(n=1600, mu=mu3, beta=beta3, pr=pr3)
fseq(finvA(q), upper = 1.21, seqlen = 26)
res = cnm(fnpvm(q), init = list(beta = 58))
plot.npvm(fnpvm(q), res$mix, res$beta, dimen = 2)

# ----



# wrapped skewed normal ---------------------------------------------------

### 1 component
n1 = 200; xi1 = pi * 0.5; omega1 = 1; alpha1 = -5; pr1 = 1
chist(rnpwsn(n1, xi1, omega1, alpha1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwsn(x, xi1, omega1, alpha1, pr1)
cdensity(f1, add = TRUE, col = 4)

# n = 100
wsn1.100.kde = sim.kde(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.100.AICc = sim.AICc(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.100.cvISE = sim.cv(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.100.cvKL = sim.cv(10, n=100, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1, cvloss="KL")

getres("wsn1", 100, mul = 10)
getres("wsn1", 100, "se", mul = 10)
getres("wsn1", 100, "numcomp")
getres("wsn1", 100, "betavalue")

# n = 400
wsn1.400.kde = sim.kde(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.400.AICc = sim.AICc(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.400.cvISE = sim.cv(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.400.cvKL = sim.cv(10, n=400, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1, cvloss="KL")

getres("wsn1", 400, mul = 10)
getres("wsn1", 400, "se", mul = 10)
getres("wsn1", 400, "numcomp")
getres("wsn1", 400, "betavalue")

# n = 1600
wsn1.1600.kde = sim.kde(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.1600.AICc = sim.AICc(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.1600.cvISE = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.1600.cvKL = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1, cvloss="KL")

getres("wsn1", 1600, mul = 10)
getres("wsn1", 1600, "se", mul = 10)
getres("wsn1", 1600, "numcomp")
getres("wsn1", 1600, "betavalue")


# test
n1 = 200; xi1 = pi * 0.5; omega1 = 0.7; alpha1 = -7; pr1 = 1
chist(rnpwsn(n1, xi1, omega1, alpha1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwsn(x, xi1, omega1, alpha1, pr1)
cdensity(f1, add = TRUE, col = 4)

wsn1.2.kde = sim.kde(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.2.AICc = sim.AICc(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.2.cvISE = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1)
wsn1.2.cvKL = sim.cv(10, n=1600, simfam="npwsn", mu=xi1, beta=omega1, alpha=alpha1, cvloss="KL")

getres("wsn1", 2, mul = 10)
getres("wsn1", 2, "se", mul = 10)
getres("wsn1", 2, "numcomp")



n1 = 200; xi1 = pi * 0.5; omega1 = 0.7; alpha1 = -7; pr1 = 1
chist(rnpwsn(n1, xi1, omega1, alpha1, pr1), nlabels = 0, nbins = 100)
f1 = function(x) dmixwsn(x, xi1, omega1, alpha1, pr1)
cdensity(f1, add = TRUE, col = 4)

# ----






