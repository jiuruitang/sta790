################################################
########## RL + Bayesian version of LM #########
########## with imputation for missing #########
################################################
cor_index = commandArgs(trailingOnly = TRUE)

if (length(cor_index) != 1)
  stop("Provide only one index")

cor_index = as.integer(cor_index)

stopifnot(!is.na(cor_index))

library(mvtnorm)
# prior
beta0 <- beta1 <- c(3:1,3:1,1,1)
alpha0 <- alpha1 <- c(rep(1,3),rep(1,3),rep(1,1),rep(1,1))

mu <- 0
tau <- 1

# Function from Mauricio
# function computes agreement levels for binary comparisons, used inside for
AgrLevBinComp <- function(x){
  same <- (x[cellinds[,1]]==x[cellinds[,2]])
  AgrLev <- 1*same
  AgrLev[!same] <- 2
  AgrLev <- as.factor(AgrLev)
  return(AgrLev)
}

# function computes agreement levels for Levenshtein comparisons, used inside for
AgrLevLevenshtein <- function(x,breaks=c(-Inf,.001,.25,.5,Inf)){
  x <- as.character(x)
  LevenshteinSim <- 1 - levenshteinSim(x[cellinds[,1]], x[cellinds[,2]])
  AgrLev <- cut(LevenshteinSim,breaks=breaks,labels=seq_len(length(breaks)-1))
  return(AgrLev)
}

gen_data <- function(my_rho=0.9){
  library(RecordLinkage)
  library(dplyr)
  data("RLdata500")
  rldf <- RLdata500
  rldf$id <- identity.RLdata500
  
  # extract all the duplicated id
  dup <- rldf[which(rldf$id %in% rldf$id[duplicated(rldf$id)]),]
  
  # arrange the records by id
  dup <- dup %>% arrange(id)
  
  # split into two files
  group <- factor(rep(1:2,50))
  dup_1 <- split(dup,group)[[1]]
  dup_2 <- split(dup,group)[[2]]
  
  # sample 12 records randomly from the not exactly matched records
  set.seed(123456)
  ss <- sample(50,12)
  fileA <- dup_1[ss,]
  fileB <- dup_2[ss,]
  
  # sample 3 records randomly from the exactly matched records
  junk <- rldf[-which(rldf$id %in% rldf$id[duplicated(rldf$id)]),]
  set.seed(1289)
  ss_m <- sample(400,3)
  exact <- junk[ss_m,]
  
  # merge matched records together, design Y in file A, and X in file B
  fileA <- rbind(fileA,exact)
  fileB <- rbind(fileB,exact)
  
  # make comparison vector weaker: change 6 of the exact match into the same 
  fileA[1:6,] <- fileA[1,]
  fileB[1:6,] <- fileB[1,]
  
  # X ~ N(0,1), Y ~ N(3+0.99*X,1)
  rho <- my_rho
  mu <- 0
  tau <- 1
  set.seed(33354)
  fileB$X <- rnorm(15,0,1)
  fileA$Y <- fileB$X*rho + rnorm(15,0,1) + 3
  
  # sample five non-matched records for each file
  junk2 <- junk[-ss_m,]
  set.seed(3435)
  ss_A <- sample(397,5)
  A_nmat <- junk2[ss_A,]
  
  junk3 <- junk2[-ss_A,]
  set.seed(43490)
  ss_B <- sample(392,5)
  B_nmat <- junk3[ss_B,]
  
  set.seed(4235767)
  B_nmat$X <- rnorm(5,0,1)
  set.seed(232303)
  known_X <- rnorm(5,0,1)
  A_nmat$Y <- known_X*rho + rnorm(5,0,1) + 3
  
  # merge together
  fileA <- rbind(fileA,A_nmat)
  fileB <- rbind(fileB,B_nmat)
  return(list(fileA,fileB))
}

fileA = gen_data(cor_index)[[1]]
fileB = gen_data(cor_index)[[2]]

# shuffle file A and file B
# fileA <- fileA[sample(dim(fileA)[1]),]
# fileB <- fileB[sample(dim(fileB)[1]),]
# idA <- fileA$id
# idB <- fileB$id

# delete ID and X,Y
record_A <- fileA[,-c(6,8,9)]
record_B <- fileB[,-c(6,8,9)]
records <- rbind(record_A,record_B)
records$by <- ifelse(records$by < 1974,0,1)
records$bd <- ifelse(records$bd<16,0,1)

n <- dim(records)[1]
n1 <- 20
n2 <- 20
cellinds <- expand.grid(1:n1,(n1+1):(n1+n2))

# computing agreement levels for binary comparisons
AgrLevBd <- AgrLevBinComp(records$bd)
AgrLevBy <- AgrLevBinComp(records$by)

# computing agreement levels for comparisons based on Levenshtein
AgrLevFN <- AgrLevLevenshtein(records$fname_c1)
AgrLevLN <- AgrLevLevenshtein(records$lname_c1)

compdata <- data.frame(i=cellinds[,1],j=cellinds[,2],GName=AgrLevFN,FName=AgrLevLN,Age=AgrLevBd,Occup=AgrLevBy)

compdata$i <- NULL
compdata$j <- NULL

nAgrLevels <- sapply(compdata,FUN=nlevels)
F <- length(nAgrLevels) # number of criteria used to compared the records
lstlevfield <- cumsum(nAgrLevels)
fstlevfield <- c(1,lstlevfield[-F]+1)
fieldInd <- rep(seq_len(F),nAgrLevels)
fieldIndParam <- fieldInd[-lstlevfield]
nBetaPriors <- sum(nAgrLevels) - F

plb1 <- rep(0,8) # prior lower bounds for binary comparisons

# next four lines compute the binary indicators from the agreement levels
xfun <- function(f) paste("(compdata[,",f,"]==",seq_len(nAgrLevels[f]),")")
expr1 <- paste(sapply(sapply(seq_len(F),xfun),FUN=paste,collapse=","),collapse=",")
expr2 <- paste("cbind(",expr1,")")
BinAgrLevels <- eval(parse(text=expr2))

BinAgrLevels[is.na(BinAgrLevels)] <- FALSE # replacing NAs by FALSE or zeroes is justified by ignorability and CI assumptions (see paper)
nBinAgrLevels <- dim(BinAgrLevels)[2]
sumsBinAgrLevels <- colSums(BinAgrLevels)

# supportive functions
which_l <- function(arr, len_arr, howmany_1s){
  which = rep(NA, (howmany_1s+1))
  j = 1
  for (i in 1:len_arr){
    if(arr[i]!=0){
      if (j < howmany_1s+1){
        which[j] = i
      }
      j = j+1
    }
  }
  return(which)
}

extract_from_to_double <- function(arr, from, to){
  len_entries = to - from + 1
  arr_entries = rep(NA,len_entries)
  for (i in 1:len_entries){
    arr_entries[i] = arr[from+i]
  }
  return(arr_entries)
}

extract_subarray_double <- function(arr, entries, len_entries){
  arr_entries = rep(NA,len_entries)
  for (i in 1:len_entries){
    arr_entries[i] = arr[entries[i]]
  }
  return(arr_entries)
}

RandomPick <- function(n,p){
  nm1 = n-1
  perm = 1:n
  # sort the prob into descending order
  temp_df <- as.data.frame(cbind(p,perm))
  temp_df <- temp_df %>% arrange(desc(p))
  p = temp_df$p
  perm = temp_df$perm
  
  # compute cumulative probabilities
  p = cumsum(p)
  
  rU = runif(1,0,p[nm1+1])
  for (i in 1:nm1){
    if (rU <= p[i])
      break;
  }
  ans = perm[i]
  return(ans)
}

# calculate likelihood in regression part
calc_llh_lm <- function(n,b,s,index){
  prob <- rep(NA,n)
  for (q in 1:n){
    prob_up <- pnorm(fileA$Y[q],mean = b[1]+fileB$X[index]*b[2],sd = s)
    prob_down <- pnorm(fileA$Y[q],mean = b[1]+b[2]*mu,sd = sqrt(s^2+b[2]^2*tau^2))
    prob[q] <- prob_up/prob_down
  }
  return(prob)
}


# Z0 doesn't link any records
# Znew <- n1:(n2+n1-1)
# freelinks <- rep(1,n1)
# Znew <- (n1+1):(n1+n2)
# Ztrue <- c(8,4,20,24,19,3,27,1,15,2,18,32,11,34,17,5,24,38,13,6)
Znew <- 1:n1
freelinks <- as.numeric(Znew > 20)

dyn.load("rtbetaC.so") # has to be in working directory

rtruncbeta <- function(alpha=2,beta=500,a=.5,iter=100){
  .C("rtbetaC", alpha = as.double(alpha), beta = as.double(beta),
     a = as.double(a), iter = as.integer(iter), x = as.double(0) )$x
}

iter <- 5000

# store the results
assign(paste("Z_",cor_index,"_imp",sep=""),matrix(NA,length(Znew),iter))
# Z_09_imp <- matrix(NA,length(Znew),iter)
assign(paste("u_",cor_index,"_imp",sep=""),matrix(NA,nBetaPriors,iter))
# u_09_imp <- matrix(NA,nBetaPriors,iter)
assign(paste("m_",cor_index,"_imp",sep=""),matrix(NA,nBetaPriors,iter))
# m_09_imp <- matrix(NA,nBetaPriors,iter)
assign(paste("betas_",cor_index,"_imp",sep=""),matrix(NA,2,iter))
# betas_09_imp <- matrix(NA,2,iter)
assign(paste("sigma_",cor_index,"_imp",sep=""),rep(NA,iter))
# sigmas_09_imp <- rep(NA,iter)

for (i in 1:iter){
  print(i)
  recinfile2withlinks <- which(Znew <= n1)
  recinfile1withlinks <- Znew[Znew <= n1]
  linked <- recinfile1withlinks + n1*(recinfile2withlinks-1)
  
  A1 <- simplify2array( lapply(seq_len(nBinAgrLevels), FUN=function(fl) sum(BinAgrLevels[linked,fl]) ) )
  A0 <- sumsBinAgrLevels - A1 
  
  revcumsumA1 <- unlist( lapply(split(A1,fieldInd), FUN=function(x) rev(cumsum(rev(x)))) )
  revcumsumA0 <- unlist( lapply(split(A0,fieldInd), FUN=function(x) rev(cumsum(rev(x)))) )
  
  A1 <- A1[-lstlevfield]
  A0 <- A0[-lstlevfield]
  
  revcumsumA1 <- revcumsumA1[-fstlevfield]
  revcumsumA0 <- revcumsumA0[-fstlevfield]
  
  ma <- A1 + alpha1
  mb <- revcumsumA1 + beta1
  ua <- A0 + alpha0
  ub <- revcumsumA0 + beta0
  
  mnewArguments <- as.matrix(cbind(ma,mb,plb1))
  mnewArgumentsList <- split(mnewArguments,seq_len(nBetaPriors))
  mnew <- simplify2array( lapply(mnewArgumentsList,FUN=function(x) rtruncbeta(alpha=x[1],beta=x[2],a=x[3],iter = 100) ) )
  unew <- rbeta(n=nBetaPriors, shape1=ua, shape2=ub)
  
  mstarnew <- unlist( lapply( split(mnew,fieldIndParam), FUN = function(x){ c(log(x),0) + c(0,cumsum(log(1-x))) } ) )
  ustarnew <- unlist( lapply( split(unew,fieldIndParam), FUN = function(x){ c(log(x),0) + c(0,cumsum(log(1-x))) } ) )
  
  mminusu <- mstarnew - ustarnew
  
  Lambda <- BinAgrLevels %*% mminusu
  
  # update regression parameters
  Y_match <- fileA$Y[Znew[Znew <= n1]]
  Y_nonmatch <- fileA$Y[Znew>n1]
  X_match <- fileB$X[Znew <= n1]
  m1 <- lm(Y_match ~ X_match)
  beta <- coefficients(m1)
  sigma_hat <- sigma(m1) 
  SSE <- sum(m1$residuals^2)
  
  # NG prior
  nu0 <- 1
  SS0 <- 1
  phi0 <- matrix(c(1,0,0,1),nrow=2,byrow = T)
  b0 <- matrix(c(3,1),nrow = 2)
  
  # NG posterior
  X <- model.matrix(m1)
  phi_n <- t(X)%*%X+phi0
  bn <- solve(phi_n)%*%(t(X)%*%X%*%beta + phi0%*%b0)
  
  # sample phi from Gamma
  gamma1 <-  length(Y_match)+nu0
  gamma2 <- SSE+SS0+t(beta)%*%t(X)%*%X%*%beta + t(b0)%*%phi0%*%b0-t(bn)%*%phi_n%*%bn
  phi <- rgamma(1,gamma1/2,gamma2/2)
  # sample beta
  new_beta_raw <- rmvnorm(1,mean = bn,sigma = solve(phi*phi_n))
  new_sigma_raw <- sqrt(1/phi)
  
  # imputation
  x_mean <- mu+new_beta_raw[2]*tau^2/(new_sigma_raw^2+new_beta_raw[2]*tau^2)*(Y_nonmatch - new_beta_raw[1]-new_beta_raw[2]*mu)
  x_var <- tau^2 - new_beta_raw[2]^2*tau^4/(new_sigma_raw^2+new_beta_raw[2]^2*tau^2)
  x_impute <- rnorm(n=x_mean, mean = x_mean,sd=sqrt(x_var))
  
  X_comb <- c(X_match,x_impute)
  Y_comb <- c(Y_match,Y_nonmatch)
  
  # Bayesian regression on more data
  m2 <- lm(Y_comb ~ X_comb)
  beta <- coefficients(m2)
  sigma_hat <- sigma(m2) 
  SSE <- sum(m2$residuals^2)
  
  # NG posterior
  X <- model.matrix(m2)
  phi_n <- t(X)%*%X+phi0
  bn <- solve(phi_n)%*%(t(X)%*%X%*%beta + phi0%*%b0)
  
  # sample phi from Gamma
  gamma1 <-  length(Y_comb)+nu0
  gamma2 <- SSE+SS0+t(beta)%*%t(X)%*%X%*%beta + t(b0)%*%phi0%*%b0-t(bn)%*%phi_n%*%bn
  phi <- rgamma(1,gamma1/2,gamma2/2)
  # sample beta
  new_beta <- rmvnorm(1,mean = bn,sigma = solve(phi*phi_n))
  new_sigma <- sqrt(1/phi)
  
  
  # update indicator Z
  nfreelinks = sum(freelinks)
  n2=as.integer(length(Znew))
  n1=as.integer(length(freelinks))
  from  = 0
  for (j in 1:n2){  #release link occupied by j, if any
    if (Znew[j] < n1){
      freelinks[Znew[j]] = 1
      # nfreelinks = nfreelinks + 1
      nfreelinks = sum(freelinks)
    }
    
    to = from + n1 -1
    
    if (nfreelinks == 0){
      Znew[j] = n1+j
    }
    else{ 
      # determine the positions of the free links
      # create a vector of length nfreelinks+1
      # which_freelinks <- which_l(freelinks,n1,nfreelinks) 
      which_freelinks <- which(freelinks == 1)
      
      # loglikelihood ratios for j and all the records from file 1
      Lambda_i <- extract_from_to_double(Lambda,from,to)
      
      # loglikelihood ratios for j and all the records from file 1 that are free to be linked
      # Lambda_i_free <- extract_subarray_double(Lambda_i,which_freelinks,nfreelinks)
      Lambda_i_free <- Lambda_i[which_freelinks]
      
      # likelihood of f(X,Y)/f(X)f(Y); q from 1 to n1; j from 1 to n2, no loop for j
      lm_llh <- calc_llh_lm(n = n2,b = new_beta,s = new_sigma,index = j)
      
      # likelihodd for j and all the records from file 1 that are free to be linked
      lm_llh_free <- lm_llh[which_freelinks]
      
      # sample new link
      Lambda_i_free = exp(Lambda_i_free)
      # add the option of non-link
      no_link = nfreelinks*(n2-n1+nfreelinks-1+1)/(n1-nfreelinks+1)
      Lambda_i_free = c(Lambda_i_free*lm_llh_free,no_link)
      nvalid_labels = nfreelinks + 1
      # which_freelinks[nfreelinks+1] = n1+j
      which_freelinks = c(which_freelinks,n1+j)
      # randomly pick one with prob Lambda_i_free
      picked_position = RandomPick(nvalid_labels, Lambda_i_free)
      Znew[j] = which_freelinks[picked_position]
      
      if (Znew[j] < n1){
        # update if j occupies any link
        freelinks[Znew[j]] = 0
        nfreelinks = nfreelinks - 1;
      }
    }
    from = from + n1
  }
  
  assign(paste("Z_",cor_index,"_imp[,",i,"]",sep=""), Znew) 
  assign(paste("m_",cor_index,"_imp[,",i,"]",sep=""), mnew)
  assign(paste("u_",cor_index,"_imp[,",i,"]",sep=""), unew)
  assign(paste("betas_",cor_index,"_imp[,",i,"]",sep=""), new_beta)
  assign(paste("sigma_",cor_index,"_imp[,",i,"]",sep=""), new_sigma)
}
save(list = c(paste("Z_",cor_index,"_imp",sep=""), 
     paste("m_",cor_index,"_imp",sep=""),
     paste("u_",cor_index,"_imp",sep=""),
     paste("betas_",cor_index,"_imp",sep=""),
     paste("sigma_",cor_index,"_imp",sep="")),
     file = paste("RL20_",cor_index,".RData",sep=""))
