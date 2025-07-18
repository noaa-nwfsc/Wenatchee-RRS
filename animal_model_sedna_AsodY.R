

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

priorR3 <- list(
  G = list(G1 = list(V = diag(2), nu = 1.002), G2 = list(V = diag(2), nu = 1.002), G3 = list(V = 1, nu = 0.002)),
  R = list(R1 = list(V = diag(2), nu = 1.002), R2 = list(V = diag(2), nu = 1.002))
)

# test version that is very short to see if it will run
model.AsodY <- MCMCglmm(Fork.off ~ as.factor(Age.off) + as.factor(Year.off) ,
                        random = ~ idh(Sex_Final.off):animal + idh(Origin_ped.off):animal + dam,
                        rcov = ~ idh(Sex_Final.off):units + idh(Origin_ped.off):units,
                        family = "gaussian",
                        pedigree = pedf, data = dfil,
                        nitt = 500000, burnin = 1000, thin = 500,
                        prior = priorR3, verbose = TRUE)
save(model.AsodY,file="model.AsodY.A_as_factor_Y_as_factor.Rdata")

summary(model.AsodY$Sol)
summary(model.AsodY$VCV)


