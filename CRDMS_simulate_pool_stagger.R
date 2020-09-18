#Huggins Closed Robust Design Multistate model simulations
#With pooling of secondary occasions
#With egress-only staggered arrival/departure during secondary occasions

#Adapated from Kery and Schaub BPA and https://sites.google.com/site/wild8390/software/simulate

#SIMULATE DATA UNDER J-S SUPERPOPULATION MODEL with Robust Design sampling

rm(list=ls())

require(RMark)


####################
# Set Parameters
N.init=100  #initial population size
n.primary=8  #number of primary sampling occassions
n.secondary=6 #number of secondary sampling occassions (per year)
S=rep(0.9,n.primary-1)    #survival probabilities for each year
p=rep(0.7,n.primary*n.secondary)      #capture probability for each secondary occasion
f=rep(0.3,n.primary-1)       #recruitment probability for each year (not estimated in CRDMS model)

n.states<-3 #number of states (including dead state)
#psiIO<-rep(0.3,n.primary-1) #transition probability for observable state (in study area) to unobservable for each year
#psiOI<-rep(0.6,n.primary-1)  #transition probability for unobservable state (not in study area) to observable
psiIO<-c(0.8,0.6,0.4,0.2,0.4,0.6,0.8) #transition probability for observable state (in study area) to unobservable for each year
psiOI<-c(0.2,0.4,0.6,0.8,0.6,0.4,0.2)  #transition probability for unobservable state (not in study area) to observable

p2.ent<-c(0.1,0.2,0.3,0.4,0,0) #arrival (and re-entry) probabilities for each secondary occasion; Assumes all arrive by 4th occasion
d<-c(0.1,0.2,0.3,0.7,0.9) #departure probabilities for each secondary interval
pers<-1-d #persistance probabilities for each secondary interval

J<-20 #number of simulations

####################
#Function to pool occasions (1-4, 5-6 within each primary occasion) in capture history
poolch<-function(CH, n.primary, n.secondary, p1=1:4, p2=5:6) {
  out <- matrix(-99, ncol=2*n.primary, nrow=dim(CH)[1]) #to store condensed capture history
  for (p in 1:n.primary) { #for each primary occasion
    out[,((p-1)*2)+1] <- apply(CH[,((p-1)*n.secondary)+p1], 1, max) #max=1 if seen at least once
    out[,((p-1)*2)+2] <- apply(CH[,((p-1)*n.secondary)+p2], 1, max)
  } #p
  return(out)
}

####################
#Function to convert capture history matrix to strings for RMark
pasty<-function(x) {
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n) {
    out[i]<-paste(x[i,],collapse="")
  }
  return(out)
}

# Function to simulate capture-recapture data under the JS Robust Design model
simul.robust.js <- function(S, p, f, psiIO, psiOI, N.init, n.primary, n.secondary, p2.ent, pers) {
  
  #get recruitment values for simulations
  #initial N
  N<-array(dim=n.primary)
  
  #B needs to be integer. get p.ent first...
  B_<-array(dim=n.primary-1)
  N[1]<-B_[1]<-N.init
  for(t in 1:(n.primary-1)) {
    B_[t+1]<-N[t]*f[t]             #number entering population at each time
    N[t+1]<-N[t]*S[t]+B_[t+1] #total population at each time, accounting for deaths and recruits
  }
  
  N.super<-round(sum(B_),0) #superpopulation (all individuals ever alive)
  p.ent<-B_/N.super        #entry probabilities
  
  B <- rmultinom(1, N.super, p.ent) #simulate number entering population at each time
  
  # 1. Survival and State transition matrix
  PSI.STATE <- array(NA, dim=c(n.states, n.states, N.super, n.primary-1))
  for (i in 1:N.super){ #for each individual...
    for (t in 1:(n.primary-1)){ #for each occasion interval...
      PSI.STATE[,,i,t] <- matrix(c(
        S[t]*(1-psiIO[t]), S[t]*psiIO[t],     1-S[t], #if observable: stays in, leaves, or dies 
        S[t]*psiOI[t],     S[t]*(1-psiOI[t]), 1-S[t], #if unobservable: arrives, stays out, or dies 
        0,           0,           1    ) #if dead: must stay dead
        , nrow=n.states, byrow=TRUE)
    } #t
  } #i
  PSI.STATE[,,i,t] #view transition matrix for last individual at last time-step
  
  # 2.Observation process matrix
  n.occasions<-n.primary*n.secondary
  PSI.OBS <- array(NA, dim=c(n.states, n.states-1, N.super, n.occasions))
  for (i in 1:N.super){ #for each individual...
    for (t in 1:n.occasions){ #for each secondary occasion...
      PSI.OBS[,,i,t] <- matrix(c(
        p[t], 1-p[t],   #if observable 
        0, 1,     #if unobservable
        0, 1    ) #if dead
        , nrow=n.states, byrow=TRUE)
    } #t
  } #i
  PSI.OBS[,,i,t] #view observation matrix for last individual at last occasion
  
  #Empty matrices to store state and captures
  CH.sur <- matrix(0, ncol=n.primary, nrow=N.super)
  CH <- CH.sur.aug <- matrix(0, ncol=n.occasions, nrow=N.super)
  
  # Define a vector with the occasion of entering the population
  ent.occ <- numeric()
  for (t in 1:n.primary){
    ent.occ <- c(ent.occ, rep(t, B[t]))
  }
  
  # Simulating survival and state transitions
  #states: 1=observable, 2=unobservable, 3=dead, 0=not yet recruited
  for (i in 1:N.super){
    CH.sur[i, ent.occ[i]] <- 1   # assume state=1 when ind. enters the pop.
    if (ent.occ[i] == n.primary) next # skip if enters pop. on last occasion
    for (t in (ent.occ[i]+1):n.primary){
      state <- which(rmultinom(1, 1, PSI.STATE[CH.sur[i,t-1],,i,t-1])==1)
      CH.sur[i,t] <- state #CH.sur holds true states of individuals at each occasion
    } #t
  } #i
  
  #fill in state info for secondary occasions
  for(i in 1:N.super) {
    x<-numeric(0)
    for(j in 1:n.primary) {
      x<-append(x,rep(CH.sur[i,j],n.secondary))
    }
    CH.sur.aug[i,]<-x
  }
  
  # Simulating arrivals and depatures during secondary occasions
  #Multiple arrivals/depatures per secondary occasion permitted
  CH.sur.arr <- CH.sur.aug #to store presence in study area
  for (i in 1:N.super) { #for each individual
    for (t in 1:n.primary){
      if (CH.sur[i,t] != 1) next # skip if not in observable state
      CH.sur.arr[i,(((t-1)*n.secondary)+1) : (((t-1)*n.secondary)+n.secondary)] <- 0 #if observable, reset secondary occasions
      ent2.occ <- which(rmultinom(1, 1, p2.ent)[,1] == 1) #simulate when ind. arrived
      CH.sur.arr[i, ((t-1)*n.secondary) + ent2.occ] <- 1   #state=1 when ind. arrived
      if (ent2.occ == n.secondary) next # skip if arrives on last occasion
      for (s in ent2.occ:n.secondary){
        ifelse (s==n.secondary, break, #break if evaluating last occasion  
                if (CH.sur.arr[i, ((t-1)*n.secondary)+s]==1){
                  CH.sur.arr[i, ((t-1)*n.secondary)+s+1] <- rbinom(1,1,pers[s]) } #if in area, persists with prob. pers
                else {
                  CH.sur.arr[i, ((t-1)*n.secondary)+s+1] <- rbinom(1,1,p2.ent[s+1]) } #if not in area, can re-enter with prob. p2.ent
        )    
      } #s
    } #t
  } #i
  
  # Simulating captures
  for (i in 1:N.super){
    for (t in 1:n.occasions){
      if (CH.sur.arr[i,t] == 0) next # skip (leave as 0) if not yet in pop., or before arrival/after departure
      event <- which(rmultinom(1, 1, PSI.OBS[CH.sur.arr[i,t],,i,t])==1)
      CH[i,t] <- event #CH holds observations, where observation of dead state may be same as unobservable state
    } #t
  } #i
  
  # Replace undetected events (event=2) with 0
  CH[CH==2] <- 0
  
  # Remove individuals never captured
  cap.sum <- rowSums(CH)
  never <- which(cap.sum == 0)
  CH <- CH[-never,]
  
  return(CH)
} #simul.robust.js



####################
# Execute simulations

#data frame to store results (egress model with pooling)
res.rdms.pool.m0<-data.frame(sim=0,estimate=-999,se=-999,lcl=-999,ucl=-999,fixed="",note="")
#data frame to store results (closed model without pooling)
res.rdms.m0<-data.frame(sim=0,estimate=-999,se=-999,lcl=-999,ucl=-999,fixed="",note="")

for (j in 1:J){ #for each simulation...
  
  #execute simulation to generate capture history
  CH<-simul.robust.js(S, p, f, psiIO, psiOI,  N.init, n.primary, n.secondary, p2.ent, pers)
  
  #######
    # Egress model with pooling
    # Pool secondary occasions
    CHp=poolch(CH, n.primary, n.secondary, p1=1:4, p2=5:6)
    
    ########
    # RMark
    
    #convert capture history to strings for RMark
    Y<-data.frame(ch=pasty(CHp))
    # Use letters to denote states
    Y$ch<-gsub("1", "A", Y$ch)
    
    #set time intervals, revised for pooling
    secondary<-c(rep(0,2-1),1)
    secondary.last<-rep(0,2-1)
    t.int<-c(rep(secondary,n.primary-1),secondary.last)
    
    # process the data for Huggins Multistate Robust Design model
    ms.data<-process.data(Y,model="HCRDMS",time.interval=t.int, strata.labels=c("A","X"))
    #define Psi parameters that are obtained by subtraction
    ms.ddl<-make.design.data(ms.data,parameters=list(Psi=list(subtract.stratum=c("A","X"))))
    #Add field to design data to fix p=0 for unobservable state
    ms.ddl$p$fix=NA
    ms.ddl$p$fix[ms.ddl$p$stratum=="X"]=0
    ms.ddl$c$fix=NA
    ms.ddl$c$fix[ms.ddl$c$stratum=="X"]=0
    
    #models
    S.dot=list(formula=~1)
    #p.dot=list(formula=~1,share=TRUE) #p=c
    #p.time=list(formula=~time,share=TRUE) #p=c varies by secondary occasion, but constant across primaries
    p.time.session=list(formula=~-1+session:time,share=TRUE) #p=c varies by time and session
    #Psi.state=list(formula=~-1+stratum:tostratum)
    Psi.state.time=list(formula=~-1+stratum:tostratum:time)
    
    #Run models
    m1=mark(ms.data,ms.ddl,model.parameters=list(S=S.dot,p=p.time.session,Psi=Psi.state.time),mlogit0=TRUE)
    unique(m1$results$real)
    
    #Store results
    out<-unique(m1$results$real);out$sim=j
    res.rdms.pool.m0<-rbind(res.rdms.pool.m0,out)
  
  
  ########
  # Closed model without pooling
  ########
  # RMark
  
  #convert capture history to strings for RMark
  Y<-data.frame(ch=pasty(CH))
  # Use letters to denote states
  Y$ch<-gsub("1", "A", Y$ch)
  
  secondary<-c(rep(0,n.secondary-1),1)
  secondary.last<-rep(0,n.secondary-1)
  t.int<-c(rep(secondary,n.primary-1),secondary.last)
  
  # process the data for Huggins Multistate Robust Design model
  ms.data<-process.data(Y,model="HCRDMS",time.interval=t.int, strata.labels=c("A","X"))
  #define Psi parameters that are obtained by subtraction
  ms.ddl<-make.design.data(ms.data,parameters=list(Psi=list(subtract.stratum=c("A","X"))))
  #Add field to design data to fix p=0 for unobservable state
  ms.ddl$p$fix=NA
  ms.ddl$p$fix[ms.ddl$p$stratum=="X"]=0
  ms.ddl$c$fix=NA
  ms.ddl$c$fix[ms.ddl$c$stratum=="X"]=0
  
  #models
  S.dot=list(formula=~1)
  #p.dot=list(formula=~1,share=TRUE) #p=c
  #p.time=list(formula=~time,share=TRUE) #p=c varies by secondary occasion, but constant across primaries
  p.time.session=list(formula=~-1+session:time,share=TRUE) #p=c varies by time and session
  #Psi.state=list(formula=~-1+stratum:tostratum)
  Psi.state.time=list(formula=~-1+stratum:tostratum:time)
  
  #Run models
  m2=mark(ms.data,ms.ddl,model.parameters=list(S=S.dot,p=p.time.session,Psi=Psi.state.time),mlogit0=TRUE)
  #unique(m2$results$real)
  
  #Store results
  out<-unique(m2$results$real);out$sim=j
  res.rdms.m0<-rbind(res.rdms.m0,out)
  
} #j

res.rdms.pool.m0
res.rdms.m0






#########
# Summaries

# Egress model with pooling
pool.res <- res.rdms.pool.m0[res.rdms.pool.m0$sim>0,] #delete placeholder row
  #add new columns of parameter
  pool.res <- transform(pool.res,
                  par = factor(substr(rownames(pool.res), 1, 3)))
  #p
  ps <- pool.res[pool.res$par=="p s",]
  #add new columns of state and time
  ps <- transform(ps,
                  s = factor(substr(rownames(ps), 4, 4)),
                  time = factor(substr(rownames(ps), 12, 13)))
  means <- with(ps, tapply(estimate, list(par,time,s), mean))
  sds <- with(ps, tapply(estimate, list(par,time,s), sd))
  cvs <- sds/means
  #No coverage or RMSE for p with open data
  means[1,,]
  cvs[1,,]
  
  #Psi
  psis <- pool.res[pool.res$par=="Psi",]
  #add new columns of state and time
  psis <- transform(psis,
                  s = factor(substr(rownames(psis), 6, 6)),
                  time = factor(substr(rownames(psis), 21, 22)))
  means <- with(psis, tapply(estimate, list(par,time,s), mean))
  sds <- with(psis, tapply(estimate, list(par,time,s), sd))
  cvs <- sds/means
  #coverage
  psis$coverage <- -999
  for (i in 1:dim(psis)[1]){  #for each row
    t <- as.integer(substr(psis[i,]$time, 2, 2)) #extract time
    if (psis[i,]$s=="A"){
      if (psis[i,]$lcl<psiIO[t] & psis[i,]$ucl>psiIO[t]){ #check coverage
        psis[i,]$coverage <- 1}
      else {psis[i,]$coverage <- 0}}
    if (psis[i,]$s=="X"){
      if (psis[i,]$lcl<psiOI[t] & psis[i,]$ucl>psiOI[t]){ #check coverage
        psis[i,]$coverage <- 1}
      else {psis[i,]$coverage <- 0}}
  } #i
  coverage <- with(psis, tapply(coverage, list(par,time,s), mean))
  #RMSE
  psis$ressq <- -999
  for (i in 1:dim(psis)[1]){  #for each row
    t <- as.integer(substr(psis[i,]$time, 2, 2)) #extract time
    if (psis[i,]$s=="A"){
      psis[i,]$ressq <- (psis[i,]$estimate - psiIO[t])^2} #calculate residual squared
    if (psis[i,]$s=="X"){
      psis[i,]$ressq <- (psis[i,]$estimate - psiOI[t])^2} #calculate residual squared
  } #i
  RMSE <- sqrt(with(psis, tapply(ressq, list(par,time,s), mean)))
  means[2,,]
  cvs[2,,]
  coverage[2,,]
  RMSE[2,,1]/psiIO #psiA:X rRMSE
  RMSE[2,,2]/psiOI #psiX:A rRMSE
  
  #S
  Ss <- pool.res[pool.res$par=="S s",]
  Ss <- transform(Ss,
                    s = factor(substr(rownames(Ss), 4, 4)),
                    time = factor(substr(rownames(Ss), 15, 16)))
  means <- with(Ss, tapply(estimate, list(par,time,s), mean))
  sds <- with(Ss, tapply(estimate, list(par,time,s), sd))
  cvs <- sds/means
  #coverage
  Ss$coverage <- -999
  for (i in 1:dim(Ss)[1]){  #for each row
    t <- as.integer(substr(Ss[i,]$time, 2, 2)) #extract time
    if (Ss[i,]$lcl<S[t] & Ss[i,]$ucl>S[t]){ #check coverage
        Ss[i,]$coverage <- 1}
      else {Ss[i,]$coverage <- 0}
  } #i
  coverage <- with(Ss, tapply(coverage, list(par,time,s), mean))
  #RMSE
  Ss$ressq <- -999
  for (i in 1:dim(Ss)[1]){  #for each row
    t <- as.integer(substr(Ss[i,]$time, 2, 2)) #extract time
    Ss[i,]$ressq <- (Ss[i,]$estimate - S[t])^2 #calculate residual squared
  } #i
  RMSE <- sqrt(with(Ss, tapply(ressq, list(par,time,s), mean)))
  means[3,,]
  cvs[3,,]
  coverage[3,,]
  RMSE[3,,]/S[1] # rRMSE; S is constant across time and states
  

#####  
# Closed model without pooling
  closed.res <- res.rdms.m0[res.rdms.m0$sim>0,] #delete placeholder row
  #add new columns of parameter
  closed.res <- transform(closed.res,
                        par = factor(substr(rownames(closed.res), 1, 3)))
  #p
  ps <- closed.res[closed.res$par=="p s",]
  #add new columns of state and time
  ps <- transform(ps,
                  s = factor(substr(rownames(ps), 4, 4)),
                  time = factor(substr(rownames(ps), 12, 13)))
  means <- with(ps, tapply(estimate, list(par,time,s), mean))
  sds <- with(ps, tapply(estimate, list(par,time,s), sd))
  cvs <- sds/means
  #No coverage or RMSE for p with open data
  means[1,,]
  cvs[1,,]
  
  #Psi
  psis <- closed.res[closed.res$par=="Psi",]
  #add new columns of state and time
  psis <- transform(psis,
                    s = factor(substr(rownames(psis), 6, 6)),
                    time = factor(substr(rownames(psis), 21, 22)))
  means <- with(psis, tapply(estimate, list(par,time,s), mean))
  sds <- with(psis, tapply(estimate, list(par,time,s), sd))
  cvs <- sds/means
  #coverage
  psis$coverage <- -999
  for (i in 1:dim(psis)[1]){  #for each row
    t <- as.integer(substr(psis[i,]$time, 2, 2)) #extract time
    if (psis[i,]$s=="A"){
      if (psis[i,]$lcl<psiIO[t] & psis[i,]$ucl>psiIO[t]){ #check coverage
        psis[i,]$coverage <- 1}
      else {psis[i,]$coverage <- 0}}
    if (psis[i,]$s=="X"){
      if (psis[i,]$lcl<psiOI[t] & psis[i,]$ucl>psiOI[t]){ #check coverage
        psis[i,]$coverage <- 1}
      else {psis[i,]$coverage <- 0}}
  } #i
  coverage <- with(psis, tapply(coverage, list(par,time,s), mean))
  #RMSE
  psis$ressq <- -999
  for (i in 1:dim(psis)[1]){  #for each row
    t <- as.integer(substr(psis[i,]$time, 2, 2)) #extract time
    if (psis[i,]$s=="A"){
      psis[i,]$ressq <- (psis[i,]$estimate - psiIO[t])^2} #calculate residual squared
    if (psis[i,]$s=="X"){
      psis[i,]$ressq <- (psis[i,]$estimate - psiOI[t])^2} #calculate residual squared
  } #i
  RMSE <- sqrt(with(psis, tapply(ressq, list(par,time,s), mean)))
  means[2,,]
  cvs[2,,]
  coverage[2,,]
  RMSE[2,,1]/psiIO #psiA:X rRMSE
  RMSE[2,,2]/psiOI #psiX:A rRMSE
  
  #S
  Ss <- closed.res[closed.res$par=="S s",]
  Ss <- transform(Ss,
                  s = factor(substr(rownames(Ss), 4, 4)),
                  time = factor(substr(rownames(Ss), 15, 16)))
  means <- with(Ss, tapply(estimate, list(par,time,s), mean))
  sds <- with(Ss, tapply(estimate, list(par,time,s), sd))
  cvs <- sds/means
  #coverage
  Ss$coverage <- -999
  for (i in 1:dim(Ss)[1]){  #for each row
    t <- as.integer(substr(Ss[i,]$time, 2, 2)) #extract time
    if (Ss[i,]$lcl<S[t] & Ss[i,]$ucl>S[t]){ #check coverage
      Ss[i,]$coverage <- 1}
    else {Ss[i,]$coverage <- 0}
  } #i
  coverage <- with(Ss, tapply(coverage, list(par,time,s), mean))
  #RMSE
  Ss$ressq <- -999
  for (i in 1:dim(Ss)[1]){  #for each row
    t <- as.integer(substr(Ss[i,]$time, 2, 2)) #extract time
    Ss[i,]$ressq <- (Ss[i,]$estimate - S[t])^2 #calculate residual squared
  } #i
  RMSE <- sqrt(with(Ss, tapply(ressq, list(par,time,s), mean)))
  means[3,,]
  cvs[3,,]
  coverage[3,,]
  RMSE[3,,]/S[1] # rRMSE; S is constant across time and states
  