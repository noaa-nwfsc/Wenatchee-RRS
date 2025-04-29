library(ggplot2)

####  script to check for potential errors in pot.par and correct them

setwd("/Users/mike.ford/Documents/MikeFord/WenSpChkProp/FRANZWensp_analysis/Git-Wenatchee-RRS")

ms = read.csv("complete genotypes June 13 2012.csv",stringsAsFactors = F)
snps = read.csv("complete_SNPS_Oct_7_2024.csv",stringsAsFactors = F)

nrow(ms)
nrow(snps)

head(ms)
names(snps)
str(ms)
ms$missing = NA
names(ms)
for(i in 1:nrow(ms))
{
  ms$missing[i] = sum(ms[i,c(2:31)]==-5)
}
summary(ms$missing)
ms$pmissing = ms$missing/30
summary(ms$pmissing)
ms$proptyped = 1-ms$pmissing
summary(ms$proptyped)

testp = read.csv("Microsat prop typed April 10 2017.csv",stringsAsFactors = F)
head(testp)
summary(testp$prop)
testp = merge(testp,ms,by="id")
head(testp)
ggplot(data=testp,aes(x=proptyped,y=prop)) + geom_point(position=position_jitter())

problems = subset(testp,prop<=.5 & proptyped >= .5)
nrow(problems)
problems$id

summary(ms)
summary(snps)
ggplot(data=snps,aes(x=missing)) + geom_histogram()
snps$proptyped = 1 - snps$missing/96
ggplot(data=snps,aes(x=proptyped)) + geom_histogram()
head(ms)
pot.par.ms = subset(ms,proptyped>=.5)
nrow(pot.par.ms)
nrow(ms)
pot.par.ms = data.frame(id = pot.par.ms$id, pot.par = 1)
head(pot.par.ms)

pot.par.snp = subset(snps,proptyped >= .5)
pot.par.snp = data.frame(id=pot.par.snp$id,pot.par = 1)
nrow(pot.par.snp)
nrow(snps)

ped = read.csv("pedigree_fitness_Tumwater_Grun_3.12.25.csv", stringsAsFactors = F)
names(ped)

years = subset(ped,select=c("id","Year.off"))
head(years)
pot.par.ms = merge(pot.par.ms,years,by="id")
pot.par.ms = subset(pot.par.ms,Year.off <=2006)

pot.par.snp = merge(pot.par.snp,years,by="id")
pot.par.snp = subset(pot.par.snp,Year.off > 2006)
pot.par = rbind(pot.par.ms,pot.par.snp)
nrow(pot.par)
nrow(pot.par.ms)
unique(duplicated(pot.par$id))
pot.par = subset(pot.par,select = c("id","pot.par"))
head(pot.par)
ped = merge(ped,pot.par,by = "id",all.x = T)
table(ped$potpar.off,ped$pot.par,useNA = 'ifany')
table(ped$pot.par,ped$noff.off,useNA = 'ifany')

names(ped)
names(ped)[145] = "pot.par.off"
ped = merge(ped,pot.par,by.x = "sire",by.y = "id",all.x = T)
names(ped)[146] = "pot.par.sire"
ped = merge(ped,pot.par,by.x = "dam",by.y = "id",all.x = T)
names(ped)[147] = "pot.par.dam"
table(ped$pot.par.off,ped$noff.off,useNA = 'ifany')
prob.off = subset(ped,noff.off >0  & is.na(pot.par.off))
nrow(prob.off)
prob.off$id
prob.off$Year.off

prob.sire = subset(ped,noff.sire >0  & is.na(pot.par.sire))
nrow(prob.sire)
prob.sire$sire
prob.sire$Year.sire

prob.dam = subset(ped,noff.dam > 0 & is.na(pot.par.dam))
nrow(prob.dam)
prob.dam$dam
prob.dam$Year.dam

ped$pot.par.off = ifelse(ped$noff.off>0,1,ped$pot.par.off) # correct the mistakes
ped$pot.par.dam = ifelse(ped$noff.dam>0,1,ped$pot.par.dam)
ped$pot.par.sire = ifelse(ped$noff.sire>0,1,ped$pot.par.sire)

# delete old columns
names(ped)
pednew = subset(ped,select = -c(potpar.off,potpar.sire,potpar.dam))
names(pednew)

write.csv(pednew,"pedigree_fitness_Tumwater_Grun_4.29.25.csv",row.names = F)
