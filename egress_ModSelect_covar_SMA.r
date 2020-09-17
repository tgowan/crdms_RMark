#SEUS SMA right whales. 1994-2015. Closed robust design multistate, egress only
#Secondary periods based on date before which all detected animals were observed (available) at least once.
#Sightings as calves excluded
#Age model: juvenile=1-8, adult=9+
rm(list=ls(all=TRUE))
library(RMark)

setwd("C:/Users/tim.gowan/Documents/AgeModel/SMA")
#setwd("E:/FWC_WAH/demography/SMA")
ms=convert.inp("whales_egress_age_2015_no-calf_SMA.inp", use.comments=TRUE, group.df=data.frame(sex=rep(c("female","male","unk"),9),age=rep(as.factor(1:9),each=3)))
dim(ms)
head(ms)
#22 primary periods each with 2 secondary occasions
t.int<-c(rep(c(0,1),21),c(0))
t.int

# process the data for multi-state model and make design data
#S=non-breeder in SEUS, V=calvingF in SEUS,  X=NoSEUS/unobservable
ms.data<-process.data(ms,model="HCRDMS",time.interval=t.int,begin.time=1994,groups=c("sex","age"), initial.age=1:9, strata.labels=c("S","V","X"))
#define Psi parameters that are obtained by subtraction
ms.ddl<-make.design.data(ms.data,parameters=list(Psi=list(subtract.stratum=c("X","X","X"))))
names(ms.ddl) #using HCRDMS (Huggins) so no f0 parameter
#view Psi at time 1 for each state, noting which parameter obtained by subtraction
ms.ddl$Psi[ms.ddl$Psi$stratum=="S"&ms.ddl$Psi$time==1994,]

#Add effort **for each secondary period** as covariate
#On-effort km in SMA
df=data.frame(esession=rep(1994:2015,2),etime=rep(1:2,each=22),
              effort=c(16321, 31007, 30245, 20284, 22099, 34139, 44363, 43546, 44253, 39545, 74364, 72741, 72609, 68017, 77884, 67422, 65754,
                       44141, 54224, 47916, 41913, 42099, 15972, 4793, 8591, 21023, 23753, 24138, 14555, 4322, 3309, 3394, 1353, 0, 2546,
                       12295, 1192, 3477, 6816, 25038, 37006, 20774, 21659, 7937))
            
start<-ms.data$begin.time
ms.ddl$p$effort<-0 #create new field and set=0 for all years
for (i in 1:dim(ms.ddl$p)[1]){ #for each row of ms.ddl$p
  session <- start + (as.numeric(ms.ddl$p$session[i])-1)     #determine year
  time <- as.numeric(ms.ddl$p$time[i])     #determine secondary period
  effort <- df$effort[df$esession==session & df$etime==time] #extract effort for that year and secondary period
  ms.ddl$p$effort[i] <- effort #update effort in design data
}
names(ms.ddl$p)[names(ms.ddl$p)=='effort']<-'effort_onkm'
ms.ddl$p[1:25,]

ms.ddl$c$effort<-0 #create new field and set=0 for all years
for (i in 1:dim(ms.ddl$c)[1]){ #for each row of ms.ddl$c
  session <- start + (as.numeric(ms.ddl$c$session[i])-1)     #determine year
  effort <- df$effort[df$esession==session & df$etime==2] #extract effort for that year and the 2nd secondary period
  ms.ddl$c$effort[i] <- effort #update effort in design data
}
names(ms.ddl$c)[names(ms.ddl$c)=='effort']<-'effort_onkm'
ms.ddl$c[1:25,]


#Add field to design data to fix p=0 for unobservable state and Calving Males/Unk
ms.ddl$p$fix=NA
ms.ddl$p$fix[ms.ddl$p$stratum=="X"]=0
ms.ddl$p$fix[ms.ddl$p$sex=="male"&ms.ddl$p$stratum=="V"]=0
ms.ddl$p$fix[ms.ddl$p$sex=="unk"&ms.ddl$p$stratum=="V"]=0
ms.ddl$p[ms.ddl$p$session==1994,]
#Add field to design data to fix c=0
ms.ddl$c$fix=NA
ms.ddl$c$fix[ms.ddl$c$stratum=="X"]=0
ms.ddl$c$fix[ms.ddl$c$sex=="male"&ms.ddl$c$stratum=="V"]=0
ms.ddl$c$fix[ms.ddl$c$sex=="unk"&ms.ddl$c$stratum=="V"]=0
ms.ddl$c[ms.ddl$c$session==1994,]

#Add design covariate to set survival different for Juves vs Adult
ms.ddl$S$young=0 #create new field and set=0 for all cases
ms.ddl$S$young[ms.ddl$S$Age<9]=1 #set=1 for calves and juves
ms.ddl$S$adult=1-ms.ddl$S$young #create new field and set=1 for adults
ms.ddl$S[ms.ddl$S$Cohort==0&ms.ddl$S$time==1994,]

#Add design covariate to set Psi different for Juves vs Adult
ms.ddl$Psi$young=0 #create new field and set=0 for all cases
ms.ddl$Psi$young[ms.ddl$Psi$Age<9]=1 #set=1 for calves and juves
ms.ddl$Psi$adult=1-ms.ddl$Psi$young #create new field and set=1 for adults
ms.ddl$Psi[ms.ddl$Psi$Cohort==0&ms.ddl$Psi$time==1994,]

#Psi for females inherently different than males/unks, so add design covariate to make it explicit
ms.ddl$Psi$fem=0 #create new field and set=0 for all cases
ms.ddl$Psi$fem[ms.ddl$Psi$sex=='female']=1 #set=1 for female
ms.ddl$Psi$nonfem=1-ms.ddl$Psi$fem #create new field and set=1 for males/unks
ms.ddl$Psi[ms.ddl$Psi$Cohort==0&ms.ddl$Psi$time==1994,]


#Add field to design data to fix impossible Psi parameters to 0
ms.ddl$Psi$fix=NA
#calving F
ms.ddl$Psi$fix[ms.ddl$Psi$stratum=="V"&ms.ddl$Psi$tostratum=="V"]=0
#Males or Unknown sex to/from Calving F
ms.ddl$Psi$fix[ms.ddl$Psi$sex=="male"&ms.ddl$Psi$tostratum=="V"]=0
ms.ddl$Psi$fix[ms.ddl$Psi$sex=="unk"&ms.ddl$Psi$tostratum=="V"]=0
ms.ddl$Psi$fix[ms.ddl$Psi$sex=="male"&ms.ddl$Psi$stratum=="V"]=0
ms.ddl$Psi$fix[ms.ddl$Psi$sex=="unk"&ms.ddl$Psi$stratum=="V"]=0
#Psi(XtoX) at time 1: optional as flag to not erroneously interpret
ms.ddl$Psi$fix[ms.ddl$Psi$stratum=="X"&ms.ddl$Psi$tostratum=="S"&ms.ddl$Psi$Time==0]=0
ms.ddl$Psi$fix[ms.ddl$Psi$stratum=="X"&ms.ddl$Psi$tostratum=="V"&ms.ddl$Psi$Time==0]=0
ms.ddl$Psi[ms.ddl$Psi$Cohort==0&ms.ddl$Psi$time==1994,]
sum(!is.na(ms.ddl$Psi$fix)) #15047 rows fixed=0

##
#Add environmental covariates for state transitions
dfc=data.frame(etime=1994:2014,
		#Gulf of Maine SST in summer of interval
		GOM=c(0.119, -0.433, -1.127, -0.853, -0.519, 0.786, 0.119, -0.312, 0.391, -0.507, -1.105, -0.391, 0.097, -0.772, -0.065, -0.100, 1.029, 0.350, 2.008, 0.898, 1.041),
		#Gulf of St Lawrence SST in summer of interval
		GSL=c(-0.682, -0.463, -0.324, -0.804, 0.046, 0.579, -0.008, 0.386, -0.667, -0.209, -0.458, 0.085, 0.885, -0.155, 0.605, 0.259, 0.686, 0.192, 1.620, 0.560, 0.660),
		#NAO index (Oct-Apr) in fall of interval through current winter														
		NAO=c(0.594, -0.494, -0.093, -0.441, 0.091, 0.786, -0.163, 0.424, -0.354, 0.283, 0.024, -0.144, 0.129, 0.286, -0.031, -1.096, -0.213, 1.086, -0.499, 0.473, 0.937),
		#NAO index (Oct-Apr), 2-yr lag
		NAOlag2=c(0.52, 1.044, 0.594, -0.494, -0.093, -0.441, 0.091, 0.786, -0.163, 0.424, -0.354, 0.283, 0.024, -0.144, 0.129, 0.286, -0.031, -1.096, -0.213, 1.086, -0.499),
    #Calanus EcoMon Gulf of Maine anomaly during interval
    CAL=c(0.1653, 0.1808, 0.0278, 0.0801, -0.0843, -0.1337, 0.0056, 0.3893, 0.0399, 0.3679, 0.5372, 0.2074, 0.1278, 0.3533,
                  0.0581, -0.0609, -0.0212, -0.1231, -0.3912, -0.0243, 0.0156),
    #Calanus EcoMon Gulf of Maine anomaly averaged during interval and prior year
    CAL2avg=c(-0.0423, 0.17305, 0.1043, 0.05395, -0.0021, -0.109, -0.06405, 0.19745, 0.2146, 0.2039, 0.45255, 0.3723, 0.1676,
                       0.24055, 0.2057, -0.0014, -0.04105, -0.07215, -0.25715, -0.20775, -0.00435))

start<-ms.data$begin.time
ms.ddl$Psi$GOM<-99 #create new fields and set=99 for all years
ms.ddl$Psi$GSL<-99
ms.ddl$Psi$NAO<-99
ms.ddl$Psi$NAOlag2<-99
ms.ddl$Psi$CAL<-99
ms.ddl$Psi$CAL2avg<-99
for (i in 1:dim(ms.ddl$Psi)[1]){ #for each row of ms.ddl$Psi
  time <- start + (as.numeric(ms.ddl$Psi$time[i])-1)     #determine year
  gom <- dfc$GOM[dfc$etime==time] #extract covariates for that year
  gsl <- dfc$GSL[dfc$etime==time]
  nao <- dfc$NAO[dfc$etime==time]
  naol2 <- dfc$NAOlag2[dfc$etime==time]
  cal <- dfc$CAL[dfc$etime==time]
  cal2av <- dfc$CAL2avg[dfc$etime==time]
  ms.ddl$Psi$GOM[i] <- gom #update value in design data
  ms.ddl$Psi$GSL[i] <- gsl
  ms.ddl$Psi$NAO[i] <- nao
  ms.ddl$Psi$NAOlag2[i] <- naol2
  ms.ddl$Psi$CAL[i] <- cal
  ms.ddl$Psi$CAL2avg[i] <- cal2av
}
ms.ddl$Psi[1:21,] 


#Update Time field to allow for quadratic trend
ms.ddl$Psi$Time=ms.ddl$Psi$Time+1 #now 1st year is 1 instead of 0



###############

memory.limit(10000)

#Run an initial model S(.)p(.)Psi(s:ts) to obtain initial values
S.dot=list(formula=~1)
p.dot=list(formula=~1, share=TRUE) #p=c
Psi.state=list(formula=~-1+ (stratum:toV) + (stratum:toS)) #X:S and S:S same for all sexes
ms.start=mark(ms.data,ms.ddl,model.parameters=list(S=S.dot,p=p.dot,Psi=Psi.state),mlogit0=TRUE)
#ms.start=mark(ms.data,ms.ddl,model.parameters=list(S=S.dot,p=p.dot,Psi=Psi.state),mlogit0=TRUE,options="SIMANNEAL")
#dimnames(model.matrix( ~-1+ (stratum:toV) + (stratum:toS), ms.ddl$Psi)) [[2]]

#View real parameter estimates
print(unique(ms.start$results$real))
get.real(ms.start,"S",se=TRUE)[get.real(ms.start,"S",se=TRUE)$time==1994 & get.real(ms.start,"S",se=TRUE)$Cohort==0 & get.real(ms.start,"S",se=TRUE)$sex=="female",]
get.real(ms.start,"p",se=TRUE)[get.real(ms.start,"p",se=TRUE)$session==1995,]
get.real(ms.start,"Psi",se=TRUE)[get.real(ms.start,"Psi",se=TRUE)$Age==1 & get.real(ms.start,"Psi",se=TRUE)$Cohort==0,]

#Derived abundance estimates
# Note 3 sexes*9 initial ages(1-9) =  27 groups, but no Unk5, Unk6, Unk7, or Unk8 in data, so 23 groups
dN <- data.frame(year=rep(1994:2015, 3*23),
                 sex=rep(c(rep('F1',22),rep('M1',22),rep('U1',22),rep('F2',22),rep('M2',22),rep('U2',22),rep('F3',22),rep('M3',22),rep('U3',22),
                           rep('F4',22),rep('M4',22),rep('U4',22),rep('F5',22),rep('M5',22),rep('F6',22),rep('M6',22),rep('F7',22),rep('M7',22),
                           rep('F8',22),rep('M8',22),rep('F9',22),rep('M9',22),rep('U9',22)),3),
                 state=c(rep('S',22*23),rep('V',22*23),rep('X',22*23)),
                 N=ms.start$results$derived$'N Population Size'[,1])
n <- rep(NA, 22); for (i in 1:length(n)){n[i] <- sum(dN$N[(dN$year-1993)==i], na.rm=T)}; plot(n~c(1994:2015),type='l')



###############
# Focus on detection

#S and Psi
S.young=list(formula=~young)
Psi.1.time=list(formula=~-1+ young:stratum:toV+young:stratum:toV:time +
                                      adult:stratum:toV+adult:stratum:toV:time +
                                      young:fem:stratum:toS+young:fem:stratum:toS:time +
                                      adult:fem:stratum:toS+adult:fem:stratum:toS:time +
                                      young:nonfem:stratum:toS+young:nonfem:stratum:toS:time +
                                      adult:nonfem:stratum:toS+adult:nonfem:stratum:toS:time ) #separate y-intercept and separate time effect for each group

#Detection
  p.dot=list(formula=~1, share=TRUE) #p=c
  p.state=list(formula=~-1+stratum, share=TRUE) #varies with state (breeders different than non-breeders)
  p.effort=list(formula=~effort_onkm, share=TRUE) #varies with effort (different for each secondary period)
  p.state.effort=list(formula=~-1+stratum+effort_onkm, share=TRUE) #varies with state, additive effort effect
  p.state.xeffort=list(formula=~stratum*effort_onkm, share=TRUE) #varies with state, interactive effort effect

#   p.time=list(formula=~-1+session:time, share=TRUE) #varies across all primary and secondary periods
#   p.state.time=list(formula=~-1+stratum+(session:time), share=TRUE) #varies with state, additive time effect
#   p.state.xtime=list(formula=~-1+stratum*(session:time), share=TRUE) #varies with state, interactive time effect
  
  #Run each model
  m.pdot=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.dot,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)
  m.pstate=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.state,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)
  m.peffort=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.effort,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)
  m.pstateeffort=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.state.effort,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)
  m.pstatexeffort=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.state.xeffort,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)

#   m.ptime=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.time,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)
#   m.pstatetime=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.state.time,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)
#   m.pstatextime=mark(ms.data,ms.ddl,model.parameters=list(S=S.young,p=p.state.xtime,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)

collect.models() #collects results from all models
#
#write.csv(collect.models()$model.table, 'prelim.AICc.mig.csv')


# Focus on survival

#p and Psi
p.state.xeffort=list(formula=~stratum*effort_onkm, share=TRUE) #varies with state, interactive effort effect
Psi.1.time=list(formula=~-1+ young:stratum:toV+young:stratum:toV:time +
                                      adult:stratum:toV+adult:stratum:toV:time +
                                      young:fem:stratum:toS+young:fem:stratum:toS:time +
                                      adult:fem:stratum:toS+adult:fem:stratum:toS:time +
                                      young:nonfem:stratum:toS+young:nonfem:stratum:toS:time +
                                      adult:nonfem:stratum:toS+adult:nonfem:stratum:toS:time ) #separate y-intercept and separate time effect for each group

#Survival
S.dot=list(formula=~1)

  #Run each model
  m.Sdotx=mark(ms.data,ms.ddl,model.parameters=list(S=S.dot,p=p.state.xeffort,Psi=Psi.1.time),mlogit0=TRUE,initial=ms.start)

collect.models() #collects results from all models
#









write.csv(collect.models()$model.table, 'prelim.AICc.csv')

#############
#Plot estimates
library(ggplot2)
#m.int <- m.pstateeffort #model of interest
unique(m.int$results$real)
m.int$results$singular
m.int$results$beta

#extract 'toV' Psi estimates
psi.v <- unique(m.int$results$real[grep('toV', rownames(m.int$results$real)), ])
psi.v <- subset(psi.v, fixed != 'Fixed')
#add new columns of state, age-class, and year
psi.v <- transform(psi.v, 
                   state = factor(substr(rownames(psi.v), 6, 6)),
                   stage = factor(ifelse(as.integer(substr(rownames(psi.v), 28, 29))>8, 'Adult', 'Juvenile'), levels = c("Juvenile", "Adult")),
                   year = as.integer(substr(rownames(psi.v), nchar(rownames(psi.v)) - 3, nchar(rownames(psi.v)))))
psi.v
g.psiv <- ggplot(psi.v, aes(x=year, y=estimate, ymin=lcl, ymax=ucl)) +
  geom_pointrange(color='slateblue') + facet_grid(state ~ stage, labeller = label_both) +
  theme_bw()
g.psiv

#extract p estimates
p.real <- unique(m.int$results$real[grep('p', rownames(m.int$results$real)), ])
p.real <- subset(p.real, fixed != 'Fixed')
#add new columns of state, secondary period, and year
p.real <- transform(p.real, 
                    state = factor(substr(rownames(p.real), 4, 4)),
                    occ = factor(substr(rownames(p.real), 22, 22)),
                    year = factor(substr(rownames(p.real), 16, 19)) )
p.real
g.preal <- ggplot(p.real, aes(x=year, y=estimate, ymin=lcl, ymax=ucl)) +
  geom_pointrange(color = 'slateblue') + facet_grid(state ~ occ, labeller = label_both) +
  theme_bw()
g.preal
###
#p-star
p.real$est.comp <- 1-p.real$estimate
pstar <- 1-(with(p.real, tapply(est.comp, list(year,state), prod)))
plot(pstar[,2], type='b', col='red', ylim=c(0,1))
points(pstar[,1], type='b', col='black')

#extract 'toS' Psi estimates
psi.s <- unique(m.int$results$real[grep('toS', rownames(m.int$results$real)), ])
psi.s <- subset(psi.s, fixed != 'Fixed')
#add new columns of state, age-class, and year
psi.s <- transform(psi.s, 
                   state = factor(substr(rownames(psi.s), 6, 6)),
                   sex = factor(substr(rownames(psi.s), 13, 13)),
                   stage = substr(rownames(psi.s), nchar(rownames(psi.s)) - 7, nchar(rownames(psi.s)) - 6),
                   year = as.integer(substr(rownames(psi.s), nchar(rownames(psi.s)) - 3, nchar(rownames(psi.s)))))
psi.s <- transform(psi.s,
                   age = factor(ifelse(substr(psi.s$stage, 1, 1)=="a" & as.integer(substr(psi.s$stage, 2, 2))<9, 
                                       'Juvenile', 'Adult'), levels = c("Juvenile", "Adult")))
psi.s
g.psis <- ggplot(psi.s, aes(x=year, y=estimate, ymin=lcl, ymax=ucl)) +
  geom_pointrange(color = 'slateblue') + facet_grid(state ~ age + sex, labeller = label_both) +
  theme_bw()
g.psis
#############