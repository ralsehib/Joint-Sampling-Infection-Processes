
# This file contains functions to simulate the joint processes: the sampling process 
# and the infection process
 
# Last modified date: 26 September 2023

# Authors: Omar De La Cruz Cabrera and Razan Alsehibani.


##### Load the required packages
library(igraph)
library(Matrix)


##### Simulation parameters
n = 200  # Number of individuals
m = 200 # Number of time steps
J = 100 # number of repetitions 


##### SIS parameters
betaN = 0.01 #Rate of transmission for unknown
gammaN = 0.002 # Rate of recovery for unknown
betaK = 0.001 # Rate of transmission for known
gammaK = 0.02 # Rate of recovery for known

##### UK parameters
alpha = 0.0001 # Rate of contact at random
delta = 0.0001 # Rate of contact through network
rao = alpha + delta # Sampling rate

##### New parameters
zeta = c(0, 0.25, 0.5, 0.75, 1) # The fraction of the random based sampling 
eta = 1 # Number of individuals being sampled per unit time
PNK = 1 # Probabiliy of meeting new unknown individual per unit time (1 means every day we will meet eta indiviudals)


AB2<-readMM(file="./matrixData/AdjBAComplex_Power_1_.mtx")
A<-AB2*1 
AAB<- graph.adjacency(A, mode="undirected", weighted=NULL) 

DiffInfec = matrix(0,nrow=m,ncol = 5)
DiffSamp = matrix(0,nrow=m,ncol = 5)
MatInfec = matrix(0, nrow=m, ncol=J)
MatSamp = matrix(0, nrow=m, ncol=J)

for (j in 1:J){ 
	for (z in 1:5){  
		s = c(1,1,rep(0,n-2) )  # initial status (1:Infected, 0:Susceptible)
		k = c(1, rep(0,n-1) ) # initial status (1:Known, 0:Unknown)
	S = matrix(0,nrow = m,ncol = n)
	S[1,] = s
	K = matrix(0,nrow = m,ncol = n)
	K[1,] = k
		for (t in 2:m){   # loop over time
           for (i in 1:n){  # loop over nodes
      # For every day starting day 2, estimate the health status of every individual on that day according to a probability proportional
      # to the status of the individual on the previous day, the infection rate, and the health status of neighbors based on a network model
				if (s[i] == 0){
					if (k[i] == 0){   
						probOfNoInfection = exp(-betaN*sum(A[i,]*ifelse(s==1,1,0)))  # "ifelse" would be needed for the SIR model
						s[i] = sample(x = 0:1, replace = FALSE, 
							size=1,
							prob = c(probOfNoInfection,1-probOfNoInfection))
        }
		
        else{
          probOfNoInfection = exp(-betaK*sum(A[i,]*ifelse(s==1,1,0)))  # "ifelse" would be needed for the SIR model
          s[i] = sample(x = 0:1, replace = FALSE, 
                        size=1,
                        prob = c(probOfNoInfection,1-probOfNoInfection))
        }
      }
      
      else{
        if (k[i] == 0){  
          probOfNoRecovery = exp(-gammaN)
          s[i] = sample(x = 0:1, replace = FALSE, 
                        size = 1,
                        prob = c(1-probOfNoRecovery,probOfNoRecovery))
        }
        else{
          probOfNoRecovery = exp(-gammaK)
          s[i] = sample(x = 0:1, replace = FALSE, 
                        size = 1,
                        prob = c(1-probOfNoRecovery,probOfNoRecovery))
        }
      }
    }  # i loop
   
 # Update the probability weight vectors for sampling a new individual 
    pnet = (!k)*(1-exp(-rao*(A %*% k)))
    pnet = pnet / sum(pnet)  ## normalized (equally likely for Unknown individuals to be sampled through network)    
    pround = !k
    pround = pround / sum(pround)  ## normalized (equally likely for Unknown individuals to be sampled randomly)
    P = zeta[z]*pround + (1-zeta[z])*pnet
    P = as.list(P)
 # The following lines are only needed when we assume that we can sample 0 or more people every day (in case we assume that we meet exactly 1 person every day, then these lines are not needed)
 # fagNewK represents the assumption either we identify (or do not identify) individual(s) on that day    
  flagNewK = sample(x = 0:1, replace = FALSE, 
                      size = 1,
                      prob = c(1-PNK,PNK))
    if (sum(k) < (n-eta) && flagNewK == 1) {
      indivSampl = sample(x = 1:n,
                          size=eta, replace = FALSE, 
                          prob = P)
      for (v in 1:eta ){
        tmp = indivSampl[v]
        k[tmp] = 1  # change the status of individual i to become known (identified)  
      }
    }  
    S[t,] = s
    K[t,] = k  
  }  # t loop 

  library(Matrix)
 # calculate the average parameter values
  DiffInfec[,z] =  rowSums(S)
  DiffSamp[,z] = rowSums(K)
  MatInfec[,j] = DiffInfec[,z]
  MatSamp[,j] = DiffSamp[,z]
 } # zeta loop 
  AvgInfec = rowMeans(MatInfec)
  AvgSamp = rowMeans(MatSamp)
} # J loop 

for ( ii in 1:5){
  color = c("black", "blue", "green", "red", "darkslategray1")
  if (ii == 1 ){
    plot(DiffInfec[,ii], type = "l", main = "Different Zeta for Two Barabási–Albert Models ", xlab = "Days", ylab = "Number of Individuals",xlim = c(0,205), ylim = c(0,n), col= color[ii])
    lines(DiffSamp[,ii], lty=2, col= color[ii]) 
  }
  
  else {
    lines(DiffInfec[,ii], col = color[ii])
    lines(DiffSamp[,ii], col = color[ii], lty = 2)  
  }
  legend("bottomright", legend = c("Zeta = 0", "Zeta = 0.25", "Zeta = 0.5", "Zeta = 0.75", "Zeta = 1"), col = color, lty=1)
}

