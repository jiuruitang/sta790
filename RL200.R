rawcor = commandArgs(trailingOnly = TRUE)

cor_index = as.numeric(rawcor)

if (length(cor_index) != 1)
  stop("Provide only one index")

cor_index = paste("0.",cor_index,sep="")

cor_index = as.numeric(cor_index)

stopifnot(!is.na(cor_index))

library(mvtnorm)
# prior
beta0 <- beta1 <- c(3:1,3:1,1,1)
alpha0 <- alpha1 <- c(rep(1,3),rep(1,3),rep(1,1),rep(1,1))

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

# Parameters of distribution of X, which assumed to be known to us
mu <- 0
tau <- 1

gen_data <- function(my_rho,N,n1,n2,n_w,my_mu=0,my_tau=1){
  # N = # of not exactly matched pairs
  # n1 = # of exactly matched pairs
  # n2 = # of non-match
  # n_w = # of weakened pairs
  library(RecordLinkage)
  library(dplyr)
  data("RLdata10000")
  rldf <- RLdata10000
  rldf$id <- identity.RLdata10000
  
  # extract all the duplicated id
  dup <- rldf[which(rldf$id %in% rldf$id[duplicated(rldf$id)]),]
  
  # arrange the records by id
  dup <- dup %>% arrange(id)
  
  # split into two files
  group <- factor(rep(1:2,1000))
  dup_1 <- split(dup,group)[[1]]
  dup_2 <- split(dup,group)[[2]]
  
  # sample 120 records randomly from the not exactly matched records
  # set.seed(123456)
  set.seed(1234)
  ss <- sample(1000,N)
  fileA <- dup_1[ss,]
  fileB <- dup_2[ss,]
  
  # sample 30 records randomly from the exactly matched records
  junk <- rldf[-which(rldf$id %in% rldf$id[duplicated(rldf$id)]),]
  #set.seed(1289)
  set.seed(3434)
  ss_m <- sample(dim(junk)[1],n1)
  exact <- junk[ss_m,]
  
  # merge matched records together, design Y in file A, and X in file B
  fileA <- rbind(fileA,exact)
  fileB <- rbind(fileB,exact)
  
  # make comparison vector weaker: change 50 of the exact match into the same 
  fileA[1:n_w,] <- fileA[1,]
  fileB[1:n_w,] <- fileB[1,]
  
  # X ~ N(0,1), Y ~ N(3+0.99*X,1)
  rho <- my_rho
  mu <- my_mu
  tau <- my_tau
  #set.seed(33354)
  set.seed(3209)
  fileB$X <- rnorm((N+n1),0,tau)
  fileA$Y <- fileB$X*rho + rnorm((N+n1),0,1) + 3
  
  # sample 50 non-matched records for each file
  junk2 <- junk[-ss_m,]
  #set.seed(3435)
  set.seed(3454693)
  ss_A <- sample(dim(junk2)[1],n2)
  A_nmat <- junk2[ss_A,]
  
  junk3 <- junk2[-ss_A,]
  #set.seed(43490)
  set.seed(43293)
  ss_B <- sample(dim(junk3)[1],n2)
  B_nmat <- junk3[ss_B,]
  
  #set.seed(4235767)
  set.seed(9894)
  B_nmat$X <- rnorm(n2,0,1)
  #set.seed(232303)
  set.seed(12950)
  known_X <- rnorm(n2,0,1)
  A_nmat$Y <- known_X*rho + rnorm(n2,0,1) + 3
  
  
  # merge together
  fileA <- rbind(fileA,A_nmat)
  fileB <- rbind(fileB,B_nmat)
  return(list(fileA,fileB))
}

AB_data = gen_data(my_rho=cor_index,N=120,n1=30,n2=50,n_w=50)
fileA = AB_data[[1]]
fileB = AB_data[[2]]

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
n1 <- 200
n2 <- 200
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

########################################
# With Imputation
########################################
Znew <- 1:n1
freelinks <- as.numeric(Znew > 200)

dyn.load("rtbetaC.so") # has to be in working directory

rtruncbeta <- function(alpha=2,beta=500,a=.5,iter=100){
  .C("rtbetaC", alpha = as.double(alpha), beta = as.double(beta),
     a = as.double(a), iter = as.integer(iter), x = as.double(0) )$x
}

iter <- 100

# store the results
assign(paste("Z_",rawcor,"_imp",sep=""),matrix(NA,length(Znew),iter))
assign(paste("u_",rawcor,"_imp",sep=""),matrix(NA,nBetaPriors,iter))
assign(paste("m_",rawcor,"_imp",sep=""),matrix(NA,nBetaPriors,iter))
assign(paste("betas_",rawcor,"_imp",sep=""),matrix(NA,2,iter))
assign(paste("sigma_",rawcor,"_imp",sep=""),rep(NA,iter))


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
  
  assign(paste("Z_",rawcor,"_imp[,",i,"]",sep=""), Znew) 
  assign(paste("m_",rawcor,"_imp[,",i,"]",sep=""), mnew)
  assign(paste("u_",rawcor,"_imp[,",i,"]",sep=""), unew)
  assign(paste("betas_",rawcor,"_imp[,",i,"]",sep=""), new_beta)
  assign(paste("sigma_",rawcor,"_imp[,",i,"]",sep=""), new_sigma)
}


#########################################
# Bayesian LM with NG Prior
# Without Imputation
#########################################
Znew <- 1:n1
freelinks <- as.numeric(Znew > n1)

# iter <- 5000


# store the results
assign(paste("Z_",rawcor,sep=""),matrix(NA,length(Znew),iter))
assign(paste("u_",rawcor,sep=""),matrix(NA,nBetaPriors,iter))
assign(paste("m_",rawcor,sep=""),matrix(NA,nBetaPriors,iter))
assign(paste("betas_",rawcor,sep=""),matrix(NA,2,iter))
assign(paste("sigma_",rawcor,sep=""),rep(NA,iter))


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
  
  assign(paste("Z_",rawcor,"[,",i,"]",sep=""), Znew) 
  assign(paste("m_",rawcor,"[,",i,"]",sep=""), mnew)
  assign(paste("u_",rawcor,"[,",i,"]",sep=""), unew)
  assign(paste("betas_",rawcor,"[,",i,"]",sep=""), new_beta)
  assign(paste("sigma_",rawcor,"[,",i,"]",sep=""), new_sigma)
}

save(list = c(paste("Z_",rawcor,"_imp",sep=""), 
              paste("m_",rawcor,"_imp",sep=""),
              paste("u_",rawcor,"_imp",sep=""),
              paste("betas_",rawcor,"_imp",sep=""),
              paste("sigma_",rawcor,"_imp",sep=""),
              paste("Z_",rawcor,sep=""), 
              paste("m_",rawcor,sep=""),
              paste("u_",rawcor,sep=""),
              paste("betas_",rawcor,sep=""),
              paste("sigma_",rawcor,sep="")),
     file = paste("RL200_",gsub("[.]","",as.character(cor_index)),".RData",sep=""))