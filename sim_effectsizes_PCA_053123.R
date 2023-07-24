#April 2023, by Doc Edge
#Goal: Simulate unlinked loci that contribute to a trait, via 
#correlated direct and indirect effects. We simulate members of two populations
#simulated by a small amount of drift (optionally zero). A parental generation is created, and
#then parental pairs are formed who each have two children.
#In the offspring generation, we conduct a sib GWAS by comparing phenotype and genotype
#within sibships, and a standard GWAS by drawing one member of each sibship
#and adjusting for a specified number of PCs.
#In a second sample, we compute polygenic scores, compute PCs, and project 
#sib and standard polygenic scores on the PCs.



library(MASS) #for drawing correlated direct and indirect effects
#from a multivariate normal




###########################################################################
###########################################################################
#Helper functions for simulation
#generate ancestral allele frequences for n.loci loci under neutral sfs with Ne = n.anc.
gen.neut.sfs <- function(n.loci, n.anc){
  probs.uw <- 1/(1:(n.anc-1))
  probs <- probs.uw/sum(probs.uw)
  cdf <- cumsum(probs)
  percs <- runif(n.loci, 0 , 1)
  get.cdfind <- function(perc){
    (sum(cdf <= perc) + 1)/n.anc 
  }
  sapply(percs, get.cdfind)
}

#drift the ancestral allele frequencies via a truncated normal.
#supply a vector of allele frequncies and a drift parameter dp, and one realization
#of the post-drift allele frequencies is provided.
drift.afs <- function(afs, dp){
  drift <- rnorm(length(afs), 0, sqrt(dp*afs*(1-afs)))
  drifted <- afs + drift
  drifted[drifted < 0] <- 0
  drifted[drifted > 1] <- 1
  return(drifted)
}


#Given two parent genotypes, draw an offspring genotype
draw.offspring.geno <- function(gens.par1, gens.par2){
  allele1 <- matrix(rbinom(length(gens.par1), 1, gens.par1/2), nrow = nrow(gens.par1))
  allele2 <- matrix(rbinom(length(gens.par2), 1, gens.par2/2), nrow = nrow(gens.par2))
  allele1 + allele2
}


#Run a GWAS for one locus.
#Gens and phens are vectors of genotypes and phenotypes.
#optionally, the user can add covariates to the regression,
#set an intercept of 0 (for sib regression), or get SEs via bootstrap
#Bootstrap is done w/ 500 bootstrap samples, and the SE is *not* used for
#calculating the p value at this moment, just for the SE estimate.
extract.lm.coeffs <- function(gens, phens, covs = NULL, int0 = FALSE, boot = FALSE){
  #if monomorphic, skip	
  if(var(gens) == 0){
    return(c(0,NA,NA,1))	
  }
  #Otherwise run regression, either without covs or with them
  if(int0 == FALSE){
    if(is.null(covs)){
      return(summary(lm(phens ~ gens))$coefficients[2,])
    }
    return(summary(lm(phens ~ gens + covs))$coefficients[2,])
  }
  if(int0 == TRUE){
    if(is.null(covs)){
      if(boot == TRUE){
        coefs <- summary(lm(phens ~ gens - 1))$coefficients[1,]
        boots <- numeric(500)
        for(i in 1:500){
          samp <- sample(1:length(gens), replace = TRUE)
          boots[i] <- sum(phens[samp]*gens[samp])/sum(gens[samp]^2)
        }
        coefs[2] <- sd(boots)
        return(coefs)
      }
      return(summary(lm(phens ~ gens - 1))$coefficients[1,])
    }
    return(summary(lm(phens ~ gens + covs - 1))$coefficients[1,])
  }
}


#Wrapper that calls extract.lm.coeffs() to estimate a bunch of GWAS
est.effs <- function(phens, gens, covs = NULL, int0 = FALSE, boot = FALSE){
  t(apply(gens, 2, extract.lm.coeffs, phens = phens, covs = covs, int0 = int0, boot = boot))
}



####################################################################
####################################################################
#Helper functions for decomposition


#Some functions to do the decomposition
#returns a centering matrix of dimension n.
center_matrix <- function(n) {
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H
}


#Compute PRS decomp at thresholds 1, .01, 5e-8 from all effect sizes
#Removing NAs prevents a crash when >0 sib-based SEs are zero
prs.decomp <- function(eff_sizes, eigenvals, eigenvecs){ 
  eigenvals * apply(eigenvecs, 2, function(x,y){sum(x*y, na.rm = TRUE)}, y = eff_sizes)^2
}

err.decomp <- function(ses, eigenvals, eigenvecs){
  eigenvals * apply(eigenvecs^2, 2, function(x,y){sum(x*y, na.rm = TRUE)}, y = ses^2)
}

#Compute PRS decomp at thresholds 1, .01, 5e-8 from all effect sizes
prs.decomp.unsq <- function(eff_sizes, eigenvals, eigenvecs){ 
  eigenvals * apply(eigenvecs, 2, function(x,y){sum(x*y, na.rm = TRUE)}, y = eff_sizes)
}

#zeros out effect size estimates (b) where p is greater than threshold
beta.p.thresh <- function(b, p, thresh){
  y <- b
  y[p >= thresh] <- 0
  y
}

#return the average squared SE of the estimates where p is less than threshold
se2.avg.p.thresh <- function(ses, p, thresh){
  mean(ses[p < thresh]^2, na.rm = TRUE)
}

#Main function to do the decomposition:
#Compute estimates of direct effect variance, strat variance, and strat-DE covariance contributions
#to variance components along each PC.
est.components <- function(gwas.b, gwas.se, sib.b, sib.se, asc.p, thresh, eigenvecs, eigenvals){
  gwas.b.thresh <- beta.p.thresh(gwas.b, asc.p, thresh)
  sib.b.thresh <- beta.p.thresh(sib.b, asc.p, thresh)
  gwas.se.thresh <- beta.p.thresh(gwas.se, asc.p, thresh)
  sib.se.thresh <- beta.p.thresh(sib.se, asc.p, thresh)
  gwas.avg.se2 <- se2.avg.p.thresh(gwas.se, asc.p, thresh)
  sib.avg.se2 <- se2.avg.p.thresh(sib.se, asc.p, thresh)
  dc.gwas <- prs.decomp(gwas.b.thresh, eigenvals, eigenvecs)
  dc.sib <- prs.decomp(sib.b.thresh, eigenvals, eigenvecs)
  dc.diff <- prs.decomp(gwas.b.thresh - sib.b.thresh, eigenvals, eigenvecs)
  proj.gwas <- prs.decomp.unsq(gwas.b.thresh, eigenvals, eigenvecs)
  proj.sib <- prs.decomp.unsq(sib.b.thresh, eigenvals, eigenvecs)
  proj.diff <- prs.decomp.unsq(gwas.b.thresh - sib.b.thresh, eigenvals, eigenvecs)
  dc.gwas.se <- err.decomp(gwas.se.thresh, eigenvals, eigenvecs)
  dc.sib.se <- err.decomp(sib.se.thresh, eigenvals, eigenvecs)
  vc.direct <- dc.sib - dc.sib.se
  vc.strat <- dc.diff - dc.sib.se - dc.gwas.se
  vc.cov <- dc.gwas - dc.diff - dc.sib + 2*dc.sib.se
  var.vc.direct <- 2*dc.sib.se^2
  var.vc.strat <- 2*(dc.sib.se + dc.gwas.se)^2
  var.vc.cov <- 4*(dc.sib.se * dc.gwas.se) + 8*dc.sib.se^2
  list(vc.direct, vc.strat, vc.cov, dc.gwas, dc.sib, dc.diff, dc.gwas.se, dc.sib.se, gwas.avg.se2, sib.avg.se2, proj.gwas, proj.sib, proj.diff, var.vc.direct, var.vc.strat, var.vc.cov, sum(gwas.b.thresh != 0))
}



#estimate gamma using a regression with intercept zero of sib vs. gwas variance components on specified PCs
#As of now, PCs have to comprise one uninterrupted interval of PC #s, default is 100:max.
get.gamma <- function(gwas.b, gwas.se, sib.b, sib.se, asc.p, thresh, eigenvecs, eigenvals, weight = FALSE, incl.lb = 100, incl.ub = length(eigenvals)){
  
  gwas.b.thresh <- beta.p.thresh(gwas.b, asc.p, thresh)
  sib.b.thresh <- beta.p.thresh(sib.b, asc.p, thresh)
  gwas.se.thresh <- beta.p.thresh(gwas.se, asc.p, thresh)
  sib.se.thresh <- beta.p.thresh(sib.se, asc.p, thresh)
  dc.gwas <- prs.decomp(gwas.b.thresh, eigenvals, eigenvecs)
  dc.sib <- prs.decomp(sib.b.thresh, eigenvals, eigenvecs)
  dc.gwas.se <- err.decomp(gwas.se.thresh, eigenvals, eigenvecs)
  dc.sib.se <- err.decomp(sib.se.thresh, eigenvals, eigenvecs)
  vc.sib <- dc.sib - dc.sib.se
  vc.gwas <- dc.gwas - dc.gwas.se
  if(weight == FALSE){
    return(sqrt(lm(vc.sib[incl.lb:incl.ub] ~ vc.gwas[incl.lb:incl.ub] - 1)$coefficients[1]))
  }
  if(weight == TRUE){
    return(sqrt(lm(vc.sib[incl.lb:incl.ub] ~ vc.gwas[incl.lb:incl.ub] - 1, weights = eigenvals[incl.lb:incl.ub])$coefficients[1]) )
  }
}


#Large wrapper function to simulate populations and estimate locus-level effect sizes.
#Simulate measured effect sizes and PCs in a gwas sample and a PRS sample
get.effectsAndPCs <- function(n.loci, n.fams.pop0, n.fams.pop1, n.loci.trait, dir.eff.var, indir.eff.var, dir.indir.cor, drift.pop0, drift.pop1, env.diff, n.fams.pop0.prs = n.fams.pop0, n.fams.pop1.prs = n.fams.pop1, correct = TRUE, boot = FALSE, n.anc = 100, scale.af = FALSE){
  #A covariance matrix for simulating direct and indirect effects
  sigma.mat.dirindir <- matrix(c(dir.eff.var, rep(dir.indir.cor*sqrt(dir.eff.var)*sqrt(indir.eff.var),2), indir.eff.var), nrow = 2)
  
  #Generate ancestral frequencies from a neutral SFS, then let them drift
  anc.afs <- gen.neut.sfs(n.loci, n.anc)
  pop0.afs <- drift.afs(anc.afs, drift.pop0)
  pop1.afs <- drift.afs(anc.afs, drift.pop1)
  
  #assign direct and (parental) indirect effects to the loci
  effs.traitloci <- mvrnorm(n.loci.trait, c(0,0), sigma.mat.dirindir) 
  effs.nontraitloci <- matrix(rep(0, 2*(n.loci - n.loci.trait) ), ncol = 2)
  effs.loci <- rbind(effs.traitloci, effs.nontraitloci)
  
  #In population 0, draw genotypes for parents (from pop0 allele freqs),
  #then draw two offspring from each pair of parents
  gens.parent1.pop0 <- matrix(rbinom(n.fams.pop0 * n.loci, 2, rep(pop0.afs, n.fams.pop0)), ncol = n.fams.pop0)
  gens.parent2.pop0 <- matrix(rbinom(n.fams.pop0 * n.loci, 2, rep(pop0.afs, n.fams.pop0)), ncol = n.fams.pop0)
  gens.sib1.pop0 <- draw.offspring.geno(gens.parent1.pop0, gens.parent2.pop0)
  gens.sib2.pop0 <- draw.offspring.geno(gens.parent1.pop0, gens.parent2.pop0)
  
  #do the same for families on which PCA/PRS decomposition will be computed
  gens.parent1.pop0.prs <- matrix(rbinom(n.fams.pop0.prs * n.loci, 2, rep(pop0.afs, n.fams.pop0.prs)), ncol = n.fams.pop0.prs)
  gens.parent2.pop0.prs <- matrix(rbinom(n.fams.pop0.prs * n.loci, 2, rep(pop0.afs, n.fams.pop0.prs)), ncol = n.fams.pop0.prs)
  gens.os.pop0.prs <- draw.offspring.geno(gens.parent1.pop0.prs, gens.parent2.pop0.prs)
  
  #Compute the genetic component (direct and indirect) of offspring
  #traits in pop 0
  indir.trait.pop0 <- t(gens.parent1.pop0 + gens.parent2.pop0) %*% effs.loci[,2]
  dir.trait.pop0.sib1 <- t(gens.sib1.pop0) %*% effs.loci[,1]
  dir.trait.pop0.sib2 <- t(gens.sib2.pop0) %*% effs.loci[,1]
  gen.trait.pop0.sib1 <- dir.trait.pop0.sib1 + indir.trait.pop0
  gen.trait.pop0.sib2 <- dir.trait.pop0.sib2 + indir.trait.pop0
  
  #same for fams in which prs is computed
  indir.trait.pop0.prs <- t(gens.parent1.pop0.prs + gens.parent2.pop0.prs) %*% effs.loci[,2]
  dir.trait.pop0.prs <- t(gens.os.pop0.prs) %*% effs.loci[,1]
  gen.trait.pop0.prs <- dir.trait.pop0.prs + indir.trait.pop0.prs
  
  #Generate genotypes and genetic components of traits in population 1
  gens.parent1.pop1 <- matrix(rbinom(n.fams.pop1 * n.loci, 2, rep(pop1.afs, n.fams.pop1)), ncol = n.fams.pop1)
  gens.parent2.pop1 <- matrix(rbinom(n.fams.pop1 * n.loci, 2, rep(pop1.afs, n.fams.pop1)), ncol = n.fams.pop1)
  gens.sib1.pop1 <- draw.offspring.geno(gens.parent1.pop1, gens.parent2.pop1)
  gens.sib2.pop1 <- draw.offspring.geno(gens.parent1.pop1, gens.parent2.pop1)
  indir.trait.pop1 <- t(gens.parent1.pop1 + gens.parent2.pop1) %*% effs.loci[,2]
  dir.trait.pop1.sib1 <- t(gens.sib1.pop1) %*% effs.loci[,1]
  dir.trait.pop1.sib2 <- t(gens.sib2.pop1) %*% effs.loci[,1]
  gen.trait.pop1.sib1 <- dir.trait.pop1.sib1 + indir.trait.pop1
  gen.trait.pop1.sib2 <- dir.trait.pop1.sib2 + indir.trait.pop1
  
  gens.parent1.pop1.prs <- matrix(rbinom(n.fams.pop1.prs * n.loci, 2, rep(pop1.afs, n.fams.pop1.prs)), ncol = n.fams.pop1.prs)
  gens.parent2.pop1.prs <- matrix(rbinom(n.fams.pop1.prs * n.loci, 2, rep(pop1.afs, n.fams.pop1.prs)), ncol = n.fams.pop1.prs)
  gens.os.pop1.prs <- draw.offspring.geno(gens.parent1.pop1.prs, gens.parent2.pop1.prs)
  indir.trait.pop1.prs <- t(gens.parent1.pop1.prs + gens.parent2.pop1.prs) %*% effs.loci[,2]
  dir.trait.pop1.prs <- t(gens.os.pop1.prs) %*% effs.loci[,1]
  gen.trait.pop1.prs <- dir.trait.pop1.prs + indir.trait.pop1.prs
  
  
  ################################################################
  #Draw environmental noise. This is entirely unshared environment;
  #there is no shared env for siblings beyond the indirect genetic effects.
  #We set the variance of the unshared env. term to be ~equal to that of
  #the variance due to direct genetic effects.
  
  env.pop0.sib1 <- rnorm(n.fams.pop0, 0, 1)
  if(sd(dir.trait.pop0.sib1) != 0){env.pop0.sib1 <- env.pop0.sib1*sd(dir.trait.pop0.sib1)}
  env.pop0.sib2 <- rnorm(n.fams.pop0, 0, 1) 
  if(sd(dir.trait.pop0.sib2) != 0){env.pop0.sib2 <- env.pop0.sib2*sd(dir.trait.pop0.sib2)}
  env.pop0.prs <- rnorm(n.fams.pop0, 0, 1) 
  if(sd(dir.trait.pop0.prs) != 0){env.pop0.prs <- env.pop0.prs*sd(dir.trait.pop0.prs)}
  env.pop1.sib1 <- (rnorm(n.fams.pop1, 0, 1) + env.diff) 
  if(sd(dir.trait.pop1.sib1) != 0){env.pop1.sib1 <- env.pop1.sib1*sd(dir.trait.pop1.sib1)}
  env.pop1.sib2 <- (rnorm(n.fams.pop1, 0, 1) + env.diff) 
  if(sd(dir.trait.pop1.sib2) != 0){env.pop1.sib2 <- env.pop1.sib2*sd(dir.trait.pop1.sib2)}
  env.pop1.prs <- rnorm(n.fams.pop1, 0, 1) 
  if(sd(dir.trait.pop1.prs) != 0){env.pop1.prs <- env.pop1.prs*sd(dir.trait.pop1.prs)}
  
  #Compute trait values for each pair of siblings in each population
  trait.pop0.sib1 <- gen.trait.pop0.sib1 + env.pop0.sib1
  trait.pop0.sib2 <- gen.trait.pop0.sib2 + env.pop0.sib2
  trait.pop0.prs <- gen.trait.pop0.prs + env.pop0.prs
  trait.pop1.sib1 <- gen.trait.pop1.sib1 + env.pop1.sib1
  trait.pop1.sib2 <- gen.trait.pop1.sib2 + env.pop1.sib2
  trait.pop1.prs <- gen.trait.pop1.prs + env.pop1.prs
  
  #Run PCA combining populations 0 and 1, but using only one sibling
  #from each sibship. This is the PCA for correcting the GWAS
  if(scale.af == FALSE){
    pc.sol <- prcomp(t(cbind(gens.sib1.pop0, gens.sib1.pop1)))
  }
  if(scale.af == TRUE){
    mat <- t(cbind(gens.sib1.pop0, gens.sib1.pop1))
    af <- colMeans(mat)/2
    mat.norm <- t((t(mat) - 2*af)/sqrt(2*af*(1-af)) )
    pc.sol <- prcomp(mat.norm)
  }
  
  #GWAS in combined pop, using only one sib from each sibship.
  if(correct==TRUE){
    GWAS.comb <- est.effs(c(trait.pop0.sib1, trait.pop1.sib1), t(cbind(gens.sib1.pop0, gens.sib1.pop1)), pc.sol$x[,1])
  }
  if(correct == FALSE){
    GWAS.comb <- est.effs(c(trait.pop0.sib1, trait.pop1.sib1), t(cbind(gens.sib1.pop0, gens.sib1.pop1)))
  }
  
  
  #sib GWAS
  sibGWAS.uncorr.diff <- est.effs(c(trait.pop0.sib1 - trait.pop0.sib2, trait.pop1.sib1 - trait.pop1.sib2), t(cbind(gens.sib1.pop0 - gens.sib2.pop0, gens.sib1.pop1 - gens.sib2.pop1)), int0 = TRUE, boot = boot)
  #sibGWAS.err <- sibGWAS.uncorr.diff - effs.loci[,1]
  
  #Run PCA combining populations 0 and 1, now using prs set
  pc.sol.gwas <- pc.sol
  rm(pc.sol)
  
  mat <- t(cbind(gens.os.pop0.prs, gens.os.pop1.prs))
  af <- colMeans(mat)/2
  if(scale.af == FALSE){
    pc.sol.prs <- prcomp(t(cbind(gens.os.pop0.prs, gens.os.pop1.prs)))
  }
  if(scale.af == TRUE){
    mat.norm <- t((t(mat) - 2*af)/sqrt(2*af*(1-af)) )
    pc.sol.prs <- prcomp(mat.norm)
  }
  
  
  #PRSs from combined GWAS
  prs.gwas.comb <- t(cbind(gens.os.pop0.prs, gens.os.pop1.prs)) %*% GWAS.comb[,1]
  
  #PRS from Sib GWAS
  prs.sib <- t(cbind(gens.os.pop0.prs, gens.os.pop1.prs)) %*% sibGWAS.uncorr.diff[,1]
  trait.prs <- c(trait.pop0.prs, trait.pop1.prs)
  
  if(scale.af == TRUE){
    GWAS.comb[,1] <- GWAS.comb[,1] * sqrt(2*af*(1-af)) 
    GWAS.comb[,2] <- GWAS.comb[,2] * sqrt(2*af*(1-af))
    sibGWAS.uncorr.diff[,1] <- sibGWAS.uncorr.diff[,1] * sqrt(2*af*(1-af)) 
    sibGWAS.uncorr.diff[,2] <- sibGWAS.uncorr.diff[,2] * sqrt(2*af*(1-af))
  }
  
  
  return(list(GWAS.comb, sibGWAS.uncorr.diff, pc.sol.prs, pc.sol.gwas, prs.gwas.comb, prs.sib, trait.prs, af))
  #returns GWAS results per locus, sib GWAS results per locus, PCA in the PRS set, the GWAS PCA (probably not needed)
  #and the PRS values for the PRS subset using both GWAS and sib effect sizes, along w/ true trait values (probably not needed)
}

###################################################################
#Run a simulation
#Do we get gamma near 1 w no indirect effects, no env. diff, 
#no structure, no ascertainment, bootstrap standard error estimates?

nsim <- 100
n.loci <- 400 #This must be at least one larger than n.loci.trait
n.fams.pop0 <- 2000 #sample sizes in each population (for both GWAS and PCA). 
n.fams.pop1 <- 2000
n.fams.pop0.prs <- 2000
n.fams.pop1.prs <- 2000
n.loci.trait <- 200 #number of loci that affect the trait
dir.eff.var <- 0.1 #variance of direct effects
indir.eff.var <-  0   #variance of indirect effects
dir.indir.cor <- 0   #correlation of direct and indirect effects
drift.pop0 <- 0 #drift time since an ancestral pop
drift.pop1 <- 0
env.diff <- 0 #difference between populations in trait 
#standard deviations due to unshared environment
thresh <- 1 #ascertainment threshold
incl.lb <- 1
incl.ub <- n.loci


# set.seed(8675309)
gammas.basic.uw <- numeric(nsim)
gammas.basic.w <- numeric(nsim)

for(i in 1:10){
  print(i)
  sim <- get.effectsAndPCs(n.loci, n.fams.pop0, n.fams.pop1, n.loci.trait, dir.eff.var, indir.eff.var, dir.indir.cor, drift.pop0, drift.pop1, env.diff, boot = TRUE, n.anc = 100, scale.af = TRUE)
  print(paste0('sim done'))
  # gwas.b, gwas.se, sib.b, sib.se, asc.p, thresh, eigenvecs, eigenvals, weight = FALSE, incl.lb = 100, incl.ub = length(eigenvals)
  write.table(sim[[1]][,1],paste0('testsims/sim.',i,'.gwas.beta.txt'))
  write.table(sim[[1]][,2],paste0('testsims/sim.',i,'.gwas.se.txt'))
  write.table(sim[[1]][,4],paste0('testsims/sim.',i,'.gwas.asc.p.txt'))
  
  write.table(sim[[2]][,1],paste0('testsims/sim.',i,'.sib.beta.txt'))
  write.table(sim[[2]][,2],paste0('testsims/sim.',i,'.sib.se.txt'))
  
  write.table(sim[[3]]$rotation,paste0('testsims/sim.',i,'.eigenvectors.txt'))
  write.table(sim[[3]]$sdev^2,paste0('testsims/sim.',i,'.eigenvalues.txt'))
  
  gammas.basic.uw[i] <- get.gamma(sim[[1]][,1], sim[[1]][,2], sim[[2]][,1], sim[[2]][,2], sim[[1]][,4], thresh, sim[[3]]$rotation, sim[[3]]$sdev^2, incl.lb = incl.lb, incl.ub = incl.ub )
  gammas.basic.w[i] <- get.gamma(sim[[1]][,1], sim[[1]][,2], sim[[2]][,1], sim[[2]][,2], sim[[1]][,4], thresh, sim[[3]]$rotation, sim[[3]]$sdev^2, incl.lb = incl.lb, incl.ub = incl.ub, weight = TRUE)
  print(i)
  print(gammas.basic.uw[i])
  print(gammas.basic.w[i])
}

write.table(gammas.basic.uw,'r.unweighted.gammas.txt')
write.table(gammas.basic.w,'r.weighted.gammas.txt')

summary(gammas.basic.uw)
summary(gammas.basic.w)

#Look at how value of gamma changes as the number of PCs included changes
#(just for most recent sim)
# gamma.omit <- numeric(399)
# gamma.omit.w <- numeric(399)
# for(i in 1:399){
#   gamma.omit[i] <-  get.gamma(sim[[1]][,1], sim[[1]][,2], sim[[2]][,1], sim[[2]][,2], sim[[1]][,4], thresh, sim[[3]]$rotation, sim[[3]]$sdev^2, incl.lb = i)
#   gamma.omit.w[i] <- get.gamma(sim[[1]][,1], sim[[1]][,2], sim[[2]][,1], sim[[2]][,2], sim[[1]][,4], thresh, sim[[3]]$rotation, sim[[3]]$sdev^2, incl.lb = i, weight = TRUE)
# }
# 
# plot(gamma.omit)
# plot(gamma.omit.w)
# 
# 
# #look at relationship between p(1-p) (x axis) and loading on first two PCs (y axis)
# plot(sim[[8]]*(1-sim[[8]]), abs(sim[[3]]$rotation[,1]))
# plot(sim[[8]]*(1-sim[[8]]), abs(sim[[3]]$rotation[,2]))
# 
# 
# 
# #Wrapper to run a whole sim
# wrap.sim <- function(nsim = 50, n.loci = 400, n.fams.pop0 = 2000, n.fams.pop1 = 2000, n.fams.pop0.prs = 2000, 
#                      n.fams.pop1.prs = 2000, n.loci.trait = 200, dir.eff.var = 0.1,
#                      indir.eff.var = 0, dir.indir.cor = 0, drift.pop0 = 0, drift.pop1 = 0, env.diff = 0, thresh = 1, incl.lb = 100, 
#                      incl.ub = n.loci, boot = FALSE, scale.af = FALSE, n.anc = 100){
#   gammas.uw <- numeric(nsim)
#   gammas.w = numeric(nsim)
#   for(i in 1:nsim){
#     sim <- get.effectsAndPCs(n.loci, n.fams.pop0, n.fams.pop1, n.loci.trait, dir.eff.var, indir.eff.var, dir.indir.cor, drift.pop0, drift.pop1, env.diff, boot = boot, n.anc = n.anc, scale.af = scale.af)
#     gammas.uw[i] <- get.gamma(sim[[1]][,1], sim[[1]][,2], sim[[2]][,1], sim[[2]][,2], sim[[1]][,4], thresh, sim[[3]]$rotation, sim[[3]]$sdev^2, incl.lb = incl.lb, incl.ub = incl.ub )
#     gammas.w[i] <- get.gamma(sim[[1]][,1], sim[[1]][,2], sim[[2]][,1], sim[[2]][,2], sim[[1]][,4], thresh, sim[[3]]$rotation, sim[[3]]$sdev^2, incl.lb = incl.lb, incl.ub = incl.ub, weight = TRUE)
#     print(gammas.uw[i])
#   }
#   return(cbind(gammas.uw, gammas.w))
# }
# 
# set.seed(8675309)
# basic <- wrap.sim(nsim = 100)
# boot <- wrap.sim(nsim = 100, boot = TRUE)
# common <- wrap.sim(nsim = 100, n.anc = 20)
# allPCs <- wrap.sim(nsim = 100, incl.lb = 1)
# scaleGRM <- wrap.sim(nsim = 100, scale.af = TRUE)
# indirs.0cor <- wrap.sim(nsim = 100, indir.eff.var = 0.1)
# indirs.poscor <- wrap.sim(nsim = 100, indir.eff.var = 0.1, dir.indir.cor = 0.7)
# indirs.negcor <- wrap.sim(nsim = 100, indir.eff.var = 0.1, dir.indir.cor = -0.7)
# asc.minus3 <- wrap.sim(nsim = 100, thresh = 1e-3)
# asc.minus6 <- wrap.sim(nsim = 100, thresh = 1e-6)
# struc <- wrap.sim(nsim=100, drift.pop0 = .01, drift.pop1 = .01)
# struc.envconf <- wrap.sim(nsim=100, drift.pop0 = .01, drift.pop1 = .01, env.diff = .1)
# mostlycausal <- wrap.sim(nsim=100, n.loci.trait = 399)
# fewcausal <- wrap.sim(nsim=100, n.loci.trait = 10)
# onecausal <- wrap.sim(nsim=100, n.loci.trait = 1)
# 
# 
# 
# 
# summaryse <- function(x){c(summary(x), sd(x)/sqrt(length(x)))}
# 
# apply(basic, 2, summaryse)
# apply(boot, 2, summaryse)
# apply(common, 2, summaryse)
# apply(allPCs, 2, summaryse)
# apply(scaleGRM, 2, summaryse)
# apply(indirs.0cor, 2, summaryse)
# apply(indirs.poscor, 2, summaryse)
# apply(indirs.negcor, 2, summaryse)
# apply(asc.minus3, 2, summaryse)
# apply(asc.minus6, 2, summaryse)
# apply(struc, 2, summaryse)
# apply(struc.envconf, 2, summaryse)
# apply(mostlycausal, 2, summaryse)
# apply(fewcausal, 2, summaryse)
# apply(onecausal, 2, summaryse)
