
# von Mises ---------------------------------------------------------------

# 6 * 3
fcdensity = function(beta=1, add=TRUE, lty=1, col=1,
                     xlim=c(-0.7,3.5), ylim=c(-0.8, 0.8)) {
  d = function(x) dmixvm(x, mu = 0, beta = beta)
  cdensity(d, add = add, nlabels = 0, main = "", xlim = xlim, ylim = ylim, 
           col = col, lty = lty)
}

fcdensity(27, FALSE)
fcdensity(9, lty = 5)
fcdensity(3, lty = 2)
fcdensity(1, lty = 3)
fcdensity(0, lty = 4)

legend(x = 2.6, y = 0.9, 
       legend = c(expression(paste(kappa, " = ", 27, sep = "")),
                  expression(paste(kappa, " = ", 9, sep = "")),
                  expression(paste(kappa, " = ", 3, sep = "")),
                  expression(paste(kappa, " = ", 1, sep = "")),
                  expression(paste(kappa, " = ", 0, sep = ""))), 
       lty = c(1,5,2,3,4), box.lty = 0, cex = 1.3, y.intersp = 1.8)

# ----



# mixture -----------------------------------------------------------------

# 7 * 5

## function to plot density curves
fcdensity = function(mu=pi/2, beta=1, pr=1, add=TRUE, lty=1, col=1,
                     xlim=c(-0.5,0.7), ylim=c(-0.85, 2.0), n=500) {
  d = function(x) dmixvm(x, mu = mu, beta = beta, pr = pr)
  cdensity(d, add = add, nlabels = 0, main = "", xlim = xlim, ylim = ylim, 
           col = col, lty = lty, n = n)
}

# data
b1 = 4; b2 = 20; b3 = 100; b4 = 500
n2 = 400; mu2 = pi * c(0.2, 0.5); beta2 = b2; pr2 = c(1, 2)/3
set.seed(1789); x2 = rnpvm(n2, mu2, beta2, pr2)
# tuning
r1 = cnm(as.npvm(x2), init = list(beta = b1), plot = "n")$mix
r2 = cnm(as.npvm(x2), init = list(beta = b2), plot = "n")$mix
r3 = cnm(as.npvm(x2), init = list(beta = b3), plot = "n")$mix
r4 = cnm(as.npvm(x2), init = list(beta = b4), plot = "n")$mix
# plot

chist(x2, nlabels = 0, nbins = 72, main = "", xlim=c(-0.5, 2.8), ylim=c(-0.6, 1.8))
fcdensity(beta = b4, mu = r4$pt, pr = r4$pr, col = 2, n = 2000)
fcdensity(beta = b3, mu = r3$pt, pr = r3$pr, lty = 5, col = 2)
fcdensity(beta = b2, mu = r2$pt, pr = r2$pr, lty = 4, col = 2)
fcdensity(beta = b1, mu = r1$pt, pr = r1$pr, lty = 3, col = 2)

legend(legend = c(expression(paste(kappa, " = ", 500, sep = "")),
                  expression(paste(kappa, " = ", 100, sep = "")),
                  expression(paste(kappa, " = ", 20, sep = "")),
                  expression(paste(kappa, " = ", 4, sep = ""))), 
       x = 1.9, y = 1.5, lty = c(1,5,4,3), box.lty = 0, cex = 1.5, y.intersp = 2, col = 2)

# ----



# simulation plot ----------------------------------------------------------

# 8 * 6.5
par(mar = rep(0, 4))
par(mfrow = c(2, 2))
cex0 = 1.2

## vM2
n2 = 400; mu2 = pi * c(0.2, 0.5); beta2 = 15; pr2 = c(1, 2)/3
set.seed(1789)
chist(rnpvm(n2, mu2, beta2, pr2), nlabels = 0, nbins = 72, main = "",
      xlim = c(-1.2, 2), ylim = c(-1, 1.5))
f2 = function(x) dmixvm(x, mu2, beta2, pr2)
cdensity(f2, add = TRUE)
text(0, -1, "(a) Mixture of 2 von Mises'", cex = cex0)

## wn3
n3 = 400; mu3 = pi * c(0.2, 0.7, 1.1); beta3 = 0.9; pr3 = c(3,1,2)/6
set.seed(1789)
chist(rnpwn(n3, mu3, beta3, pr3), nlabels = 0, nbins = 72, main = "",
      xlim = c(-1.2, 2), ylim = c(-1, 1.5))
f3 = function(x) dmixwn(x, mu3, beta3, pr3)
cdensity(f3, add = TRUE)
text(0, -1, "(b) Mixture of 3 wrapped normals", cex = cex0)

## wC2
n2 = 400; mu2 = pi * c(0.2, 0.5); beta2 = 0.8; pr2 = 1:2/3
fwc2 = function(x) dmixwc(x, mu2, beta2, pr2)
set.seed(1642)
chist(rnpwc(n2, mu2, beta2, pr2), nlabels = 0, nbins = 72, main = "",
      xlim = c(-1.2, 2), ylim = c(-1, 1.5))
cdensity(fwc2, add = TRUE)
text(0, -1, "(c) Mixture of 2 wrapped Cauchy's", cex = cex0)

## wsn1
n1 = 400; xi1 = pi * 0.5; omega1 = 1; alpha1 = -5
fwsn1 = function(x) dmixwsn(x, xi1, omega1, alpha1)
set.seed(1789)
chist(rnpwsn(n1, xi1, omega1, alpha1), nlabels = 0, nbins = 72, main = "",
      xlim = c(-1.2, 2), ylim = c(-1, 1.5))
cdensity(fwsn1, add = TRUE)
text(0, -1, "(d) Wrapped skewed normal", cex = cex0)

par(mfrow = c(1, 1))

# ----



# real data setup ---------------------------------------------------------

?Turtles_radians
?HurricanesGulfofMexico1971to2008
?ProteinsAAA

## CircNNTSR::turtle
data(Turtles_radians)
turtle = Turtles_radians
turtle = 2 * pi - turtle + pi / 2
turtle[turtle >= 2 * pi] = turtle[turtle >= 2 * pi] - 2 * pi
length(turtle)

## CircNNTSR::protein
data(ProteinsAAA)
dim(ProteinsAAA)
protein = ProteinsAAA$V2
length(protein)

## CircNNTSR::hurricane2
data(HurricanesGulfofMexico1971to2008)
hurricane = HurricanesGulfofMexico1971to2008$V1 * 2 * pi
hurricane = 2 * pi - hurricane + pi / 2
hurricane[hurricane >= 2 * pi] = hurricane[hurricane >= 2 * pi] - 2 * pi
length(hurricane)

# ----



# real data plot  ---------------------------------------------------------

# 6 * 6
par(mar = rep(0, 4))

# turtle
rad = 0.59
chist(turtle, nbins = 72, nlabels = 0, main = "", radius = rad, 
      xlim = c(-1, 1.4), ylim = c(-0.75, 1.1))
f1 = function(x) dmixvm(x, turtle.all.AICc$mix$pt, turtle.all.AICc$beta, turtle.all.AICc$mix$pr)
f2 = function(x) dmixvm(x, turtle, turtle1.all.pi, 1/length(turtle))
cdensity(f2, radius = rad, add = TRUE, lty = 2)
cdensity(f1, radius = rad, add = TRUE)
p = 0.45; text(c(p, 0, -p, 0), c(0, -p, 0, p), c("E", "S", "W", "N"), cex = 1.8)

# protein
rad = 0.75
chist(protein, nbins = 72, nlabels = 0, main = "", radius = rad, 
      xlim = c(-0.7, 2.7), ylim = c(-1, 1))
f1 = function(x) dmixvm(x, protein.all.cvISE$mix$pt, protein.all.cvISE$beta, protein.all.cvISE$mix$pr)
f2 = function(x) dmixvm(x, protein, protein.all.pi, 1/length(protein))
cdensity(f2, radius = rad, add = TRUE, lty = 2)
cdensity(f1, radius = rad, add = TRUE, n = 2000)
p = 0.55; text(c(p, 0, -p, 0), c(0, p, 0, -p), cex = 1.5,
               labels = expression(0, frac(pi,2), pi, frac(3*pi,2)))

# hurricane
rad = 0.5
chist(hurricane, nbins = 72, nlabels = 0, main = "", radius = rad,
      xlim = c(-1.35, 0.55), ylim = c(-1.0, 0.5))
f1 = function(x) dmixvm(x, hurricane.all.cvISE$mix$pt, hurricane.all.cvISE$beta, hurricane.all.cvISE$mix$pr)
f2 = function(x) dmixvm(x, hurricane, hurricane.all.lscv, 1/length(hurricane))
cdensity(f2, radius = rad, add = TRUE, lty = 2)
cdensity(f1, radius = rad, add = TRUE)
p = 0.39; text(c(p, 0, -p, 0), c(0, -p, 0, p), c(6, 12, 18, 0), cex = 1.8)

# ----



# real data tuning ---------------------------------------------------------------

# results
turtle.all.AICc$beta
xtable(rbind(turtle.all.AICc$mix$pt, turtle.all.AICc$mix$pr), 
       align = rep("l", 3), digits = 3)

protein.all.cvISE$beta
xtable(rbind(protein.all.cvISE$mix$pt, protein.all.cvISE$mix$pr), 
       align = rep("l", 14), digits = 3)

hurricane.all.cvISE$beta
# hurricane.all.cvISE$mix$pt = hurricane.all.cvISE$mix$pt / 2 / pi * 24
xtable(rbind(hurricane.all.cvISE$mix$pt / 2 / pi * 24, hurricane.all.cvISE$mix$pr), 
       align = rep("l", 6), digits = 3)
round((hurricane.all.cvISE$mix$pt / 2 / pi * 24) %% 1 * 60)

## CircNNTSR::turtle
getres.real("turtle")
getres.real("turtle", "se")

turtle.all.AICc = tune("mde", "AICc", turtle, seqlen=40, upper=1.05)
# turtle.all.cvISE = tune("mde", "cv", turtle, seqlen=40, upper=1.05)
# turtle.all.cvKL = tune("mde", "cv", turtle, seqlen=40, upper=1.05, cvloss="KL")

# turtle.all.lcv = tune("kde", "LCV", turtle)
# turtle.all.lscv = tune("kde", "LSCV", turtle)
turtle.all.pi = tune("kde", "pi", turtle)
# turtle.all.rt = tune("kde", "rt", turtle)


## CircNNTSR::protein
getres.real("protein")
getres.real("protein", "se")

# protein.all.AICc = tune("mde", "AICc", protein, seqlen=100, upper=1.045)
protein.all.cvISE = tune("mde", "cv", protein, seqlen=100, upper=1.045)
# protein.all.cvKL = tune("mde", "cv", protein, seqlen=100, upper=1.045, cvloss="KL")

# protein.all.lcv = tune("kde", "LCV", protein)
# protein.all.lscv = tune("kde", "LSCV", protein)
protein.all.pi = tune("kde", "pi", protein)
# protein.all.rt = tune("kde", "rt", protein)


## CircNNTSR::hurricane2
getres.real("hurricane")
getres.real("hurricane", "se")

# hurricane.all.AICc = tune("mde", "AICc", hurricane, seqlen=50, upper=1.05)
hurricane.all.cvISE = tune("mde", "cv", hurricane, seqlen=50, upper=1.05)
# hurricane.all.cvKL = tune("mde", "cv", hurricane, seqlen=50, upper=1.05, cvloss="KL")

# hurricane.all.lcv = tune("kde", "LCV", hurricane)
hurricane.all.lscv = tune("kde", "LSCV", hurricane)
# hurricane.all.pi = tune("kde", "pi", hurricane)
# hurricane.all.rt = tune("kde", "rt", hurricane)

# ----









