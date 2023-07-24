
library(MASS) #for drawing correlated direct and indirect effects
#from a multivariate normal

nsim <- 1
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


set.seed(8675309)

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

drift.afs <- function(afs, dp){
  
  set.seed(8675309)
  drift <- rnorm(length(afs), 0, sqrt(dp*afs*(1-afs)))
  drifted <- afs + drift
  drifted[drifted < 0] <- 0
  drifted[drifted > 1] <- 1
  return(drifted)
}

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

get.effectsAndPCs <- function(n.loci, n.fams.pop0, n.fams.pop1, n.loci.trait, dir.eff.var, indir.eff.var, dir.indir.cor, drift.pop0, drift.pop1, env.diff, n.fams.pop0.prs = n.fams.pop0, n.fams.pop1.prs = n.fams.pop1, correct = TRUE, boot = FALSE, n.anc = 100, scale.af = FALSE){
  #A covariance matrix for simulating direct and indirect effects
  sigma.mat.dirindir <- matrix(c(dir.eff.var, rep(dir.indir.cor*sqrt(dir.eff.var)*sqrt(indir.eff.var),2), indir.eff.var), nrow = 2)
  
  #Generate ancestral frequencies from a neutral SFS, then let them drift
  anc.afs <- gen.neut.sfs(n.loci, n.anc)
  write.table(anc.afs,'R.and.afs.txt')
  pop0.afs <- drift.afs(anc.afs, drift.pop0)
  pop1.afs <- drift.afs(anc.afs, drift.pop1)
  
  #assign direct and (parental) indirect effects to the loci
  effs.traitloci <- mvrnorm(n.loci.trait, c(0,0), sigma.mat.dirindir) 
  effs.nontraitloci <- matrix(rep(0, 2*(n.loci - n.loci.trait) ), ncol = 2)
  effs.loci <- rbind(effs.traitloci, effs.nontraitloci)
  
  #In population 0, draw genotypes for parents (from pop0 allele freqs),
  #then draw two offspring from each pair of parents
  set.seed(8675309)
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
    print(dim(pc.sol$x))
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


# set.seed(8675309)
gammas.basic.uw <- numeric(nsim)
gammas.basic.w <- numeric(nsim)

for(i in 1:10){
  sim <- get.effectsAndPCs(n.loci, n.fams.pop0, n.fams.pop1, n.loci.trait, dir.eff.var, indir.eff.var, dir.indir.cor, drift.pop0, drift.pop1, env.diff, boot = TRUE, n.anc = 100, scale.af = TRUE)
  print(paste0('sim done'))
  # gwas.b, gwas.se, sib.b, sib.se, asc.p, thresh, eigenvecs, eigenvals, weight = FALSE, incl.lb = 100, incl.ub = length(eigenvals)
  write.table(sim[[1]][,1],paste0('testsims/sim.',i,'.gwas.beta.txt'))
  write.table(sim[[1]][,2],paste0('testsims/sim.',i,'.gwas.se.txt'))
  write.table(sim[[1]][,4],paste0('testsims/sim.',sim,'.gwas.asc.p.txt'))
  
  write.table(sim[[2]][,1],paste0('testsims/sim.',sim,'.sib.beta.txt'))
  write.table(sim[[2]][,2],paste0('testsims/sim.',sim,'.sib.se.txt'))
  
  write.table(sim[[3]]$rotation,paste0('testsims/sim.',sim,'.eigenvectors.txt'))
  write.table(sim[[3]]$sdev^2,paste0('testsims/sim.',sim,'.eigenvalues.txt'))
  
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
