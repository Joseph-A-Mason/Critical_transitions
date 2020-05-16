setwd("")

##Set number of age intervals (bins)
nclass <- 20

##Fill matrix with calculated number of sand (active) or soil 
##(stable) episodes in each age interval. An episode is one
##period of sand deposition or soil formation at one individual
##site, documented by stratigraphy and age control
m <- matrix(scan('sandsoil.txt'), nrow = nclass, byrow=TRUE)

##function rep.col replicates list 
rep.col <- function(x, n){
matrix(rep(x, each=n), ncol=n, byrow=TRUE)
}


##set number of simulations, typically 10,000
nsims <- 10000


##replicate first column of matrix m, listing number of sand episodes in  
##each interval, to fill sandrep, a matrix with nsims identical columns
sandrep <- rep.col(m[1:nclass, 1], nsims)

##create a zero vector with the same dimensions as sandrep
sandshift <- mat.or.vec(nclass, nsims)

##create a matrix representing probabilities that true age of an episode of activity
##or stability at a site places it in classes up to four younger or up to four older 
##than the one it has been placed in based on the controlling OSL age(s) for that
##episode. This matrix is read from a file containing the values of a cumulative
##normal distribution at the boundaries of the four younger and four older bins.
##The normal distribution used is derived from the average of numerical ages that
##fall within that bin and their average 1-sigma errors.In effect, it will be 
##assumed that the probability that epsiode actually belongs in a younger or older
##age class equals the probability that a typical OSL age in a particular age class 
##actually represents an event in a younger or older age class
shiftprobs <- matrix(scan('shiftprobs.txt'), nrow = nclass, byrow=TRUE)

##nsims simulations of the age distribution of sand (active) intervals are
##now generated.For each age class, sand interval are either left in this age
##class or moved to older or younger classes, based on a random number and
##the normal distribution values listed in the shiftprobs matrix. "Moved" means 
##the value of sandshift for that simulation and age class being processed is 
##reduced by one,and the value of sandshift for the age class the interval moves 
##to is increased by one. In effect, sandshift records the shifts to be applied 
##to the original distribution of sand sites by age class, in each simulation.
##intervals cannot be moved below the youngest age class or above the oldest
for(i in 1:nsims) {
  for(j in 1:nclass) {
    if(sandrep[j,i] > 0) {
      nmoves <- sandrep[j,i]
      for(h in 1:nmoves) {
        trp <- runif(1)
        if(trp < shiftprobs[j,1]) {
          if(j > 4)
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j-4,i] <- sandshift[j-4,i] + 1
          } 
        }
        if(trp < shiftprobs[j,2] & trp >= shiftprobs[j,1]){
          if(j > 3) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j-3,i] <- sandshift[j-3,i] + 1
          }
        }
        if(trp < shiftprobs[j,3] & trp >= shiftprobs[j,2]){
          if(j>2) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j-2,i] <- sandshift[j-2,i] + 1
          }
        }
        if(trp < shiftprobs[j,4] & trp >= shiftprobs[j,3]){
          if(j>1) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j-1,i] <- sandshift[j-1,i] + 1
          }
        }
        if(trp > shiftprobs[j,5] & trp <= shiftprobs[j,6]){
          if(j < nclass) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j+1,i] <- sandshift[j+1,i] + 1
          } 
        }
        if(trp > shiftprobs[j,6] & trp <= shiftprobs[j,7]){
          if(j < (nclass-1)) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j+2,i] <- sandshift[j+2,i] + 1
          }
        }
        if(trp > shiftprobs[j,7] & trp <= shiftprobs[j,8]){
          if(j < (nclass-2)) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j+3,i] <- sandshift[j+3,i] + 1
          } 
        }
        if(trp > shiftprobs[j,8] & trp <= 1){
          if(j < (nclass-3)) {
            sandshift[j,i] <- sandshift[j,i] - 1
            sandshift[j+4,i] <- sandshift[j+4,i] + 1
        } 
      }
    }
  }
}

##the sandshift matrix is used to generate nsims different simulations
##by adding to or subtracting from the original numbers of sand intervals
##in each age class.
for (g in 1:nclass) {
  for (h in 1:nsims) {
      sandrep[g, h] <- sandshift[g, h] + sandrep[g,h]
  } 
}

##replicate second column of matrix m, listing number of soil episodes in  
##each interval, to fill soilrep, a matrix with nsims identical columns
soilrep <- rep.col(m[1:nclass, 2], nsims)

##create a zero vector with the same dimensions as soilrep
soilshift <-mat.or.vec(nclass, nsims)

##nsims simulations of the age distribution of soil (stable) intervals are
##now generated.For each age class, soil interval are either left in this age
##class or moved to older or younger classes, based on a random number and
##the normal distribution values listed in the shiftprobs matrix. "Moved" means 
##the value of soilshift for that simulation and age class being processed is 
##reduced by one,and the value of soilshift for the age class the interval moves 
##to is increased by one. In effect, soilshift records the shifts to be applied 
##to the original distribution of soil sites by age class, in each simulation.
##intervals cannot be moved below the youngest age class or above the oldest
for(i in 1:nsims) {
  for(j in 1:nclass) {
    if(soilrep[j,i] > 0) {
      nmoves <- soilrep[j,i]
      for(h in 1:nmoves) {
        trp <- runif(1)
        if(trp < shiftprobs[j,1]) {
          if(j > 4) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j-4,i] <- soilshift[j-4,i] + 1
          } 
        }
        if(trp < shiftprobs[j,2] & trp >= shiftprobs[j,1]){
          if(j > 3) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j-3,i] <- soilshift[j-3,i] + 1
          }
        }
        if(trp < shiftprobs[j,3] & trp >= shiftprobs[j,2]){
          if(j>2) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j-2,i] <- soilshift[j-2,i] + 1
          }
        }
        if(trp < shiftprobs[j,4] & trp >= shiftprobs[j,3]){
          if(j>1) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j-1,i] <- soilshift[j-1,i] + 1
          }
        }
        if(trp > shiftprobs[j,5] & trp <= shiftprobs[j,6]){
          if(j < nclass) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j+1,i] <- soilshift[j+1,i] + 1
          } 
        }
        if(trp > shiftprobs[j,6] & trp <= shiftprobs[j,7]){
          if(j < (nclass-1)) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j+2,i] <- soilshift[j+2,i] + 1
          }
        }
        if(trp > shiftprobs[j,7] & trp <= shiftprobs[j,8]){
          if(j < (nclass-2)) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j+3,i] <- soilshift[j+3,i] + 1
          } 
        }
        if(trp > shiftprobs[j,8] & trp <= 1){
          if(j < (nclass-3)) {
            soilshift[j,i] <- soilshift[j,i] - 1
            soilshift[j+4,i] <- soilshift[j+4,i] + 1
          } 
        }
      }
    }
  }
}

##the soilshift matrix is used to generate nsims different simulations
##by adding to or subtracting from the original numbers of soil intervals
##in each age class.
for (g in 1:nclass) {
  for (h in 1:nsims) {
      soilrep[g, h] <- soilshift[g, h] + soilrep[g,h]
  }
}

##the matrix percstable is created to store the values of percent
##stable in each age class for each simulation
percstable <- matrix(, nrow = nclass, ncol = nsims)

##the percent stable episodes out of all episodes in each age class
##and each simulation is calculated and used to fill percstable
for (g in 1:nclass) {
  for (h in 1:nsims) {
    if(soilrep[g,h] + sandrep[g,h] > 0) {
      percstable[g, h] <- soilrep[g, h]/(sandrep[g, h] + soilrep[g,h])
    } else {
      percstable[g,h] <- 0
    }
  }
}

##matrix ci is created to store statistics for nsims simulations
ci <- matrix(, nrow = nclass, ncol = 4)

##2.5 percentile, 97.5 percentile, mean, and standard deviation of the values of 
##percent stable for each age class across nsims simulations is calculated and 
##used to fill ci
for (k in 1:nclass) {
  ci[k, 1] <- quantile(percstable[k, 1:nsims], probs = c(0.025))
  ci[k, 2] <- quantile(percstable[k, 1:nsims], probs = c(0.975))
  ci[k, 3] <- quantile(percstable[k, 1:nsims], probs = c(0.50))
  ci[k, 4] <- mean(percstable[k, 1:nsims])
}

ci

