
########################################################################
#12/7/2024
#Goal: Draw figures and conduct simulations to explore the effect
#of changing case fraction on case-control GWAS, with the perspective
#of using the bounds on r^2 in terms of allele frequencies to understand
#the effects.

#Notation follows Paye & Edge, 2024.


##############################################################3
#Helper functions

#Function that takes in parameters p (allele freq of A), 
#b (disease freq for A genotype), gamma (disease freq for a genotype),
#c (degree to which cases are oversampled).
#Returns a 2x2 matrix of joint probabilities of carrying each of the two 
#genotypes and of being a case vs. control. Haploid case.
get.cont.probs.hap <- function(p, b, gamma, c){
  d <- b*p + gamma*(1-p)
  cont.a <- p*(1-b)*(1-c*d)/(1-d)
  cont.A <- (1-d-b*p)*(1-c*d)/(1-d)
  case.a <- c*b*p
  case.A <- c*(d-b*p)
  matrix(c(cont.a, case.a, cont.A, case.A), nrow = 2)
}

#Function that takes in parameters p (allele freq of A), 
#b (disease freq for AA genotype), gamma (disease freq for aa genotype),
#c (degree to which cases are oversampled), and h (dominance coefficient).
#Returns a 2x3 matrix of joint probabilities of carrying each of the three 
#genotypes and of being a case vs. control. Diploid case.
get.cont.probs.dip <- function(p, b, gamma, c, h){
  d <- b*p^2 +(h*b + (1-h)*gamma)*2*p*(1-p) + gamma*(1-p)^2
  cont.aa <- p^2 * (1-b)*(1-c*d)/(1-d)
  cont.Aa <- 2*p*(1-p)*(1-h*b-gamma*(1-h))*(1-c*d)/(1-d)
  cont.AA <- (1-p)^2*(1-gamma)*(1-c*d)/(1-d)
  case.aa <- p^2*b*c
  case.Aa <- 2*p*(1-p)*(h*b + gamma*(1-h))*c
  case.AA <- (1-p)^2*gamma*c
  matrix(c(cont.aa, case.aa, cont.Aa, case.Aa, cont.AA, case.AA),nrow = 2)
}

#Use the Brandt-Snedecor formula to compute lambda, the chisq effect size,
#from a contingency table.
get.lam <- function(cont.probs.tab){
  tab <- cont.probs.tab
  q <- rowSums(tab)[2]
  fi <- colSums(tab)
  qi <- tab[2,]/fi
  sum(fi*(qi-q)^2)/(q*(1-q))
}


#Given a contingency probability table, sample size n,
#and type I error rate alpha, return the expected power of a standard 
#pearson chisq test for association.
get.power <- function(cont.probs.tab, n, alpha){
  tab <- cont.probs.tab
  lam <- get.lam(tab)
  ncp <- n*lam
  df <- ncol(tab)-1
  crit <- qchisq(1-alpha, df, ncp=0)
  1-pchisq(crit, df, ncp)
}


#get lambda (effect size) values as a function of c for fixed p, b, and gamma
lamcurve.hap <- function(p, b, gamma){
  d <- b*p + gamma*(1-p)
  c.vals <- seq(0, 1/d, length.out = 10000)
  lams <- numeric(length(c.vals))
  for(i in 1:length(c.vals)){
    tab <- get.cont.probs.hap(p, b, gamma, c.vals[i])
    lams[i] <- get.lam(tab)
  }
  return(cbind(c.vals * d, lams))
}


#get diploid lambda (effect size) values as a function of c for 
#fixed p, b, gamma, h.
lamcurve.dip <- function(p, b, gamma, h){
  d <- b*p^2 +(h*b + (1-h)*gamma)*2*p*(1-p) + gamma*(1-p)^2
  c.vals <- seq(0, 1/d, length.out = 10000)
  lams <- numeric(length(c.vals))
  for(i in 1:length(c.vals)){
    tab <- get.cont.probs.dip(p, b, gamma, c.vals[i], h)
    lams[i] <- get.lam(tab)
  }
  return(cbind(c.vals * d, lams))
}


#Get diploid expected X2 power as a function of c for fixed p, b, gamma, h,
#and n. also requires T1 error rate
powcurve.dip <- function(p, b, gamma, h, n = 10000, alpha = 5e-8){
  d <- b*p^2 +(h*b + (1-h)*gamma)*2*p*(1-p) + gamma*(1-p)^2
  c.vals <- seq(0, 1/d, length.out = 10000)
  pow <- numeric(length(c.vals))
  for(i in 1:length(c.vals)){
    tab <- get.cont.probs.dip(p, b, gamma, c.vals[i], h)
    pow[i] <- get.power(tab, n, alpha)
  }
  return(cbind(c.vals * d, pow))
}

#pc.1 <- powcurve.dip(p = 0.1, b = .75, gamma = 0.1, h = 0, n = 1000)
#plot(pc.1, type = "l", bty = "n", xlab = "Case fraction", ylab = "power", las = 1)
#lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")


#Simulate a contingency table of data given the probabilities
#in cont.tab. Returns n observations. If fix.cases = TRUE, then
#the ratio of cases will match cont.tab as close to exactly as possible.
#Works for haploid or diploid contingency tables
sim.cont.tab <- function(cont.tab, n, fix.cases = TRUE){
  if(fix.cases == FALSE){
    counts <- rmultinom(1, n, as.numeric(cont.tab))
    return(matrix(counts, nrow = 2))
  }
  if(fix.cases == TRUE){
    n.case <- round(sum(cont.tab[2,]) * n)
    counts.cont <- as.numeric(rmultinom(1, n - n.case, cont.tab[1,]))
    counts.case <- as.numeric(rmultinom(1, n.case, cont.tab[2,]))
    return(rbind(counts.cont, counts.case))
  }
}

#gentab <- sim.cont.tab(cont.tab, 1000)
#chisq.test(gentab)
#prop.trend.test(geno.tab[2,], colSums(geno.tab))


#performs a logistic regression using a contingency table of
#genotypes by case status as input. Does a probit if the linkf argument is
#set to "probit"
glm.from.genotab <- function(geno.tab, linkf = "logit"){
  x <- numeric(0)
  y <- numeric(0)
  for(i in 1:ncol(geno.tab)){
    xi <- rep(i-1, sum(geno.tab[,i]))
    x <- c(x,xi)
    yi <- c(rep(0, geno.tab[1,i]), rep(1, geno.tab[2,i]))
    y <- c(y, yi)
  }
  glm(y ~ x, family = binomial(linkf))
}

#summary(glm.from.genotab(geno.tab, linkf = "probit"))


##Given a contingency table of probabilities and a sample size,
#returns empirical power (at level alpha) from nsim 
#simulations of chisq test
cspow.conttab <- function(cont.tab, n, nsim, alpha = 5e-8){
  ps <- numeric(nsim)
  for(i in 1:nsim){
    gentab <- sim.cont.tab(cont.tab, n)
    suppressWarnings(ps[i] <- chisq.test(gentab[,colSums(gentab) > 0])$p.value)
  }
  mean(ps < alpha, na.rm = TRUE)
}



#Given a contingency table of probabilities and a sample size,
#returns empirical power (at level alpha) from nsim 
#simulations of chisq test, armitage
#trend test, logistic regression, and probit regression.
#In haploid case, excludes armitage. 
pows.conttab <- function(cont.tab, n, nsim, alpha = 5e-8){
  if(ncol(cont.tab) == 2){
    ps <- matrix(nrow = nsim, ncol = 4)
  }
  if(ncol(cont.tab) == 3){
    ps <- matrix(nrow = nsim, ncol = 4)
  }
  for(i in 1:nsim){
    gentab <- sim.cont.tab(cont.tab, n)
    suppressWarnings(cs.p <- chisq.test(gentab)$p.value)
    logit.p <- summary(glm.from.genotab(gentab))$coefficients[2,4]
    probit.p <- summary(glm.from.genotab(gentab), linkf = "probit")$coefficients[2,4]
    if(ncol(cont.tab) == 2){
      ps[i,] <- c(cs.p, logit.p, probit.p)
    }
    if(ncol(cont.tab) == 3){
      cat.p <- prop.trend.test(gentab[2,], colSums(gentab))$p.value
      ps[i,] <- c(cs.p, cat.p, logit.p, probit.p)
    }
  }
  colMeans(ps < alpha, na.rm = TRUE)
}


#plot expected X2 power and realized X2 power as a function of c for fixed p,
#b, gamma, h, and n, and type 1 error rate alpha. computes empirical power 
#for sample fractions in fracs using nsim simulations. Draws +=2*se bars around
#each empirical power.
#if maxpow is not null, then chooses n such that the maximum theoretical
#power is equal to maxpower. In this case, the "n" parameter is overridden.
powfig.cs <- function(p, b, gamma, h, n = 10000, alpha = 5e-8, 
                      fracs = (1:9)/10, nsim = 1000, maxpow = NULL){
  if(!is.null(maxpow)){
    lc <- lamcurve.dip(p, b, gamma, h)
    lam.max <- max(lc[,2], na.rm = TRUE)
    crit.chisq <- qchisq(1-alpha, 2)
    crit.ncp <- optim(1, function(ncp){abs(pchisq(crit.chisq, 2, ncp)-(1-maxpow))}, method = "Brent", lower = 0, upper = 100)
    n <- round(crit.ncp$par/lam.max)
    print(n)
  }
  pc <- powcurve.dip(p, b, gamma, h, n, alpha)
  d <- b*p^2 +(h*b + (1-h)*gamma)*2*p*(1-p) + gamma*(1-p)^2
  cs <- fracs/d
  pows <- numeric(length(fracs))
  for(i in 1:length(fracs)){
    ct <- get.cont.probs.dip(p, b, gamma, cs[i], h)
    pows[i] <- cspow.conttab(ct, n, nsim, alpha)
  }
  plot(pc, type = "l", bty = "n", xlab = "Case fraction", ylab = "Power", las = 1, ylim = c(0,1))
  points(fracs, pows, pch = 20)
  for(i in 1:length(fracs)){
    se <- sqrt(pows*(1-pows)/nsim)
    lines(c(fracs[i], fracs[i]), c(pows[i] - 2*se[i], pows[i] + 2*se[i]))
  }
  lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
  text(paste("n =", as.character(n)), x = 0.85, y = 0.95)
}



#plot expected X2 power and realized power from all tests 
#as a function of c for fixed p,
#b, gamma, h, and n, and type 1 error rate alpha. computes empirical power 
#for sample fractions in fracs using nsim simulations. Draws +=2*se bars around
#each empirical power.
#if maxpow is not null, then chooses n such that the maximum theoretical
#power is equal to maxpower. In this case, the "n" parameter is overridden.
powfig.allmethods <- function(p, b, gamma, h, n = 10000, alpha = 5e-8, 
                      fracs = (1:9)/10, nsim = 1000, maxpow = NULL, pal = 
                        c('#a6cee3','#1f78b4','#b2df8a','#33a02c'), leg = TRUE ){
  if(!is.null(maxpow)){
    lc <- lamcurve.dip(p, b, gamma, h)
    lam.max <- max(lc[,2], na.rm = TRUE)
    crit.chisq <- qchisq(1-alpha, 2)
    crit.ncp <- optim(1, function(ncp){abs(pchisq(crit.chisq, 2, ncp)-(1-maxpow))}, method = "Brent", lower = 0, upper = 100)
    n <- round(crit.ncp$par/lam.max)
    print(n)
  }
  pc <- powcurve.dip(p, b, gamma, h, n, alpha)
  d <- b*p^2 +(h*b + (1-h)*gamma)*2*p*(1-p) + gamma*(1-p)^2
  cs <- fracs/d
  pows <- matrix(nrow = length(fracs), ncol = 4)
  for(i in 1:length(fracs)){
    ct <- get.cont.probs.dip(p, b, gamma, cs[i], h)
    pows[i,] <- pows.conttab(ct, n, nsim, alpha)
  }
  plot(pc, type = "l", bty = "n", xlab = "Case fraction", ylab = "Power", las = 1, ylim = c(0,1))
  points(fracs - .015, pows[,1], pch = 20, col = pal[1])
  points(fracs - .005, pows[,2], pch = 20, col = pal[2])
  points(fracs + .005, pows[,3], pch = 20, col = pal[3])
  points(fracs + .015, pows[,4], pch = 20, col = pal[4])
  se <- sqrt(pows*(1-pows)/nsim)
  for(i in 1:length(fracs)){
    lines(c(fracs[i] - .015, fracs[i] - .015), c(pows[i,1] - 2*se[i,1], pows[i,1] + 2*se[i,1]), col = pal[1])
    lines(c(fracs[i] - .005, fracs[i] - .005), c(pows[i,2] - 2*se[i,2], pows[i,2] + 2*se[i,2]), col = pal[2])
    lines(c(fracs[i] + .005, fracs[i] + .005), c(pows[i,3] - 2*se[i,3], pows[i,3] + 2*se[i,3]), col = pal[3])
    lines(c(fracs[i] + .015, fracs[i] + .015), c(pows[i,4] - 2*se[i,4], pows[i,4] + 2*se[i,4]), col = pal[4])
  }
  lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
  text(paste("n =", as.character(n)), x = 0.85, y = 0.95)
  if(leg == TRUE){
    legend("right", bty = "n", col = pal, pch = 20, 
           legend = c("Chi-squared", "Trend test", "Logistic", "Probit"), cex = 0.8)
  }
}

#powfig.allmethods(p = 0.1, b = 0.3, gamma = 0.05, h = 1/2, maxpow = 0.7)


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

#Figure 2: lambda as a function of the case fraction in the 
#haploid case at varying effect sizes 

lc.1 <- lamcurve.hap(p = 0.2, b=0.1, gamma = 0.05)
lc.2 <- lamcurve.hap(p = 0.2, b=0.8, gamma = 0.05)
lc.3 <- lamcurve.hap(p = 0.2, b=1, gamma = 0.05)

pdf("fig2.pdf", width = 7, height = 2.5)
par(mfrow = c(1,3), mar = c(4,4,2,2))
plot(lc.1, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("A", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.2, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("B", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.3, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("C", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
dev.off()

#######################################################################
#######################################################################
#Figure 3: lambda as a function of the case fraction in the 
#diploid case at varying effect sizes and dominance

lc.1 <- lamcurve.dip(p = 0.1, b = 0.1, gamma = 0.05, h = 0)
lc.2 <- lamcurve.dip(p = 0.1, b = 0.1, gamma = 0.05, h = 0.5)
lc.3 <- lamcurve.dip(p = 0.1, b = 0.1, gamma = 0.05, h = 1)
lc.4 <- lamcurve.dip(p = 0.1, b = 0.8, gamma = 0.05, h = 0)
lc.5 <- lamcurve.dip(p = 0.1, b = 0.8, gamma = 0.05, h = 0.5)
lc.6 <- lamcurve.dip(p = 0.1, b = 0.8, gamma = 0.05, h = 1)
lc.7 <- lamcurve.dip(p = 0.1, b = 1, gamma = 0.05, h = 0)
lc.8 <- lamcurve.dip(p = 0.1, b = 1, gamma = 0.05, h = 0.5)
lc.9 <- lamcurve.dip(p = 0.1, b = 1, gamma = 0.05, h = 1)

pdf("fig3.pdf", width = 7, height = 7)
par(mfrow = c(3,3), mar = c(4,4,2,2))
plot(lc.1, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("A", adj = 0)
mtext("h=0")
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.2, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("B", adj = 0)
mtext("h=1/2")
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.3, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("C", adj = 0)
mtext("h=1")
mtext("b=0.1", side = 4)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.4, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("D", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.5, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("E", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.6, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("F", adj = 0)
mtext("b=0.8", side = 4)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.7, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("G", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.8, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("H", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.9, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("I", adj = 0)
mtext("b=1", side = 4)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
dev.off()


#############################################################
############################################################
#Figure 4: minor risk alleles vs. minor protective alleles


lc.1 <- lamcurve.dip(p = 0.1, b = 0.0125, gamma = 0.1, h = 0)
lc.2 <- lamcurve.dip(p = 0.1, b = 0.0125, gamma = 0.1, h = 0.5)
lc.3 <- lamcurve.dip(p = 0.1, b = 0.0125, gamma = 0.1, h = 1)
lc.4 <- lamcurve.dip(p = 0.1, b = 0.8, gamma = 0.1, h = 0)
lc.5 <- lamcurve.dip(p = 0.1, b = 0.8, gamma = 0.1, h = 0.5)
lc.6 <- lamcurve.dip(p = 0.1, b = 0.8, gamma = 0.1, h = 1)


pdf("fig4.pdf", width = 7, height = 4.5)
par(mfrow = c(2,3), mar = c(4,4,2,2))
plot(lc.1, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("A", adj = 0)
mtext("h=0")
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.2, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("B", adj = 0)
mtext("h=1/2")
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.3, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("C", adj = 0)
mtext("h=1")
mtext("minor allele protective", side = 4)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.4, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("D", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.5, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("E", adj = 0)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
plot(lc.6, type = "l", bty = "n", xlab = "Case fraction", ylab = expression(paste(lambda, " (effect size)")), las = 1)
mtext("F", adj = 0)
mtext("major allele protective", side = 4)
lines(c(0.5, 0.5), c(0,1), lty = 2, col = "gray")
dev.off()


#############################################################
############################################################
#Figure 5: predicted power roughly maps to actual power

set.seed(8675309)
pdf("fig5.pdf", width = 7, height = 7)
par(mfrow = c(3,3), mar = c(4,4,2,2))
powfig.cs(p = 0.1, b = 0.05, gamma = 0.02, h = 0, maxpow = 0.9)
mtext("A", adj = 0)
mtext("h=0")
powfig.cs(p = 0.1, b = 0.05, gamma = 0.02, h = 0.5, maxpow = 0.9)
mtext("B", adj = 0)
mtext("h=1/2")
powfig.cs(p = 0.1, b = 0.05, gamma = 0.02, h = 1, maxpow = 0.9)
mtext("C", adj = 0)
mtext("h=1")
mtext("b=0.05", side = 4)
powfig.cs(p = 0.1, b = 0.2, gamma = 0.02, h = 0, maxpow = 0.9)
mtext("D", adj = 0)
powfig.cs(p = 0.1, b = 0.2, gamma = 0.02, h = 0.5, maxpow = 0.9)
mtext("E", adj = 0)
powfig.cs(p = 0.1, b = 0.2, gamma = 0.02, h = 1, maxpow = 0.9)
mtext("F", adj = 0)
mtext("b=0.2", side = 4)
powfig.cs(p = 0.1, b = .4, gamma = 0.02, h = 0, maxpow = 0.9)
mtext("G", adj = 0)
powfig.cs(p = 0.1, b = .4, gamma = 0.02, h = 0.5, maxpow = 0.9)
mtext("H", adj = 0)
powfig.cs(p = 0.1, b = .4, gamma = 0.02, h = 1, maxpow = 0.9)
mtext("I", adj = 0)
mtext("b=.4", side = 4)
dev.off()


#############################################################
############################################################
#Figure 6: predicted and empirical power for recessive alleles with large effect sizes

set.seed(8675309)
pdf("fig6.pdf", width = 7, height = 4.5)
par(mfrow = c(2,3), mar = c(4,4,2,2))
powfig.cs(p = 0.1, b = 0.8, gamma = 0.1, h = 0, maxpow = 0.9)
mtext("A", adj = 0)
mtext("b=0.8")
powfig.cs(p = 0.1, b = 0.9, gamma = 0.1, h = 0, maxpow = 0.9)
mtext("B", adj = 0)
mtext("b=0.9")
powfig.cs(p = 0.1, b = 1, gamma = 0.1, h = 0, maxpow = 0.9)
mtext("C", adj = 0)
mtext("b=1")
mtext("p=0.1", side = 4)
powfig.cs(p = 0.01, b = 0.8, gamma = 0.1, h = 0, maxpow = 0.9)
mtext("D", adj = 0)
powfig.cs(p = 0.01, b = 0.9, gamma = 0.1, h = 0, maxpow = 0.9)
mtext("E", adj = 0)
powfig.cs(p = 0.01, b = 1, gamma = 0.1, h = 0, maxpow = 0.9)
mtext("F", adj = 0)
mtext("p=.01", side = 4)
dev.off()



######################################################################
###################################################################
#Figure 7: Power for other methods
set.seed(8675309)
pdf("fig7.pdf", width = 7, height = 4.5)
par(mfrow = c(2,3), mar = c(4,4,2,2))
powfig.allmethods(p = 0.1, b = 0.1, gamma = 0.02, h = 0, maxpow = 0.7, leg = FALSE)
mtext("A", adj = 0)
mtext("h=0")
powfig.allmethods(p = 0.1, b = 0.1, gamma = 0.02, h = 0.05, maxpow = 0.7, leg = FALSE)
mtext("B", adj = 0)
mtext("h=1/20")
powfig.allmethods(p = 0.1, b = 0.1, gamma = 0.02, h = 1/2, maxpow = 0.7, leg = FALSE)
mtext("C", adj = 0)
mtext("h=1/2")
mtext("b=0.1, p=.1", side = 4)
powfig.allmethods(p = 0.01, b = 0.8, gamma = 0.02, h = 0, maxpow = 0.7, leg = TRUE)
mtext("D", adj = 0)
powfig.allmethods(p = 0.01, b = 0.8, gamma = 0.02, h = 0.05, maxpow = 0.7, leg = FALSE)
mtext("E", adj = 0)
powfig.allmethods(p = 0.01, b = 0.8, gamma = 0.02, h = 1/2, maxpow = 0.7, leg = FALSE)
mtext("F", adj = 0)
mtext("b=0.8, p=.01", side = 4)

dev.off()
