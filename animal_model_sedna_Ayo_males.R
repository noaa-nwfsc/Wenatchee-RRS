

####### script for running animal model on sedna 

ssh mford@sedna.nwfsc2.noaa.gov

screen -h 10000
screen -ls
screen -d -R
screen -d -r 4118683


### note use ctrl-a esc to enter scrolling mode and then just esc to exit scrolling mode

ls
cd /home/mford/animal_models
ls
squeue
sinfo --Node
#srun -c 20 --nodelist node33 --time 30-0 --pty /bin/bash
srun --time 20-0 --pty /bin/bash
##srun -c 24 -p himem --pty /bin/bash

module load R
R
q("no")
exit
exit
exit


library(MCMCglmm)
library(pedtricks)
setwd("/home/mford/animal_models")


#### load data

d = read.csv("pedigree_fitness_Tumwater_Grun_4.29.25.csv",stringsAsFactors = F,na.strings = "NA")
d$id = as.factor(d$id)
d$dam = as.factor(d$dam)
d$sire = as.factor(d$sire)
d = d[,c(3,1,2,4:ncol(d))]
names(d)[1] = "animal"

#### filter data

ped = d
ped = subset(ped,Origin_ped.off == "Hatchery" | Origin_ped.off == "Wild")
ped = subset(ped,!is.na(Origin_ped.off))
ped = subset(ped,!is.na(Origin_ped.dam))
ped = subset(ped,!is.na(Origin_ped.sire))
ped = subset(ped,Sex_Final.off == "Male" | Sex_Final.off == "Female")
ped = subset(ped,GRun.off == "Spring" )
ped = subset(ped,Final_Status_ped.off == "Broodstock" | Final_Status_ped.off == "Spawning Grounds")
ped = subset(ped,Year.off <=2018 & Year.off >= 2008) # after this hatchery fish were not sampled
# limit to neither parent NA
ped = subset(ped, !is.na(sire) & !is.na(dam)) ### try not doing this

##### format pedigree

pedf = fix_ped(ped)
#dfil = subset(d,animal %in% pedf$id)
dfil = ped

table(dfil$Origin_ped.off,dfil$Year.off)


##### set up and run model

##### set up and run fkl model - 
##### male only

dmal = subset(dfil,Sex_Final.off=="Male")
pvar = var(dmal$Fork.off,na.rm=T)

priorR2 <- list(
  G = list(G1 = list(V = diag(2)*pvar*.25, nu = .2), G2 = list(V = pvar*.25, nu = 0.2), G3 = list(V=pvar*.25,nu=.2)),
  R = list(R2 = list(V = diag(2)*pvar*.5, nu = .2))
)

# fork length model, age fixed effects and origin and year as random effects
model.Adyo_mal <- MCMCglmm(Fork.off ~ as.factor(Age.off) ,
                                random = ~ idh(Origin_ped.off):animal + dam + Year.off ,rcov = ~idh(Origin_ped.off):units,
                                family = "gaussian",
                                pedigree = pedf, data = dmal,
                                nitt = 500000, burnin = 10000, thin = 500,
                                prior = priorR2, verbose = TRUE)
save(model.Adyo_mal,file="model.Adyo_mal.Rdata")

##### set up and run fkl model - 
##### female only

dmal = subset(dfil,Sex_Final.off=="Female")
pvar = var(dmal$Fork.off,na.rm=T)

priorR2 <- list(
  G = list(G1 = list(V = diag(2)*pvar*.25, nu = .2), G2 = list(V = pvar*.25, nu = 0.2), G3 = list(V=pvar*.25,nu=.2)),
  R = list(R2 = list(V = diag(2)*pvar*.5, nu = .2))
)

# fork length model, age fixed effects and origin and year as random effects
model.Adyo_fem <- MCMCglmm(Fork.off ~ as.factor(Age.off) ,
                                random = ~ idh(Origin_ped.off):animal + dam + Year.off ,rcov = ~idh(Origin_ped.off):units,
                                family = "gaussian",
                                pedigree = pedf, data = dmal,
                                nitt = 500000, burnin = 10000, thin = 500,
                                prior = priorR2, verbose = TRUE)
save(model.Adyo_fem,file="model.Adyo_fem.Rdata")


