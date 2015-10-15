# Battery.R
# This script computes the optimal battery size given a power consumption profile and maximum cable capacity
# 
# Log:
# - Start  16-01-15: This script does all kinds of battery simulations (I have no idea which yet, but it has something to do with linear optimization).
# - Update 19-01-15: This script now calculates the optimal battery size to overcome the MSR peak load, isn't it awesome?!
# - Update 15-10-15: Made the script fit for a HAN workshop by breaking it and scaling it down a bit. (I have removed the most important line.)
# 
# By Werner van Westering MSc.
# Start 16-01-2015

##Load packages
print("--Loading packages--")
library(slam)
library(Rglpk) 
library(tictoc)
library(doSNOW)
library(parallel)
library(data.table)

# writeLines(c(""), "log.txt")
# cl<-makeCluster(7) 
# registerDoSNOW(cl)
# sink("log.txt")
# sink("log.txt", append=TRUE,type = "message")
# 
# A = foreach(iter = 1:100, .packages='tictoc', .combine=rbind, .verbose=FALSE) %dopar% {
#   tic()
#   sink("log.txt")
#   sink("log.txt", append=TRUE,type = "message")
#   cat(paste("Starting iteration",toc(),"\n"))
#   return(iter)
# }
# stopCluster(cl)
# sink(type = "message")
# sink()

##Load data
gc(verbose=FALSE)
drive = 'E:/'
path = paste0(drive,"1. Programmeerwerk/Bottum Up Analyse/2. Data")
setwd(paste0(path,"/7. Output"))
print("--Loading data--")
load("Connections_NH_v2.RData")
rm(baseloadperHLD,baseloadperOSLD)
# load("Batteryinput.Rdata")
gc(verbose=FALSE)

print("--Calculating reference scenario--")
##Calculate reference scenario
Assetii = 49
ii      = nscenarios
profile = AllKVtechprofiles %*% KVScenariosperMSR[Assetii,,ii] + AllGVtechprofiles %*% GVScenariosperMSR[Assetii,,ii] 
profile = profile[-((length(profile)-96-1):length(profile))]
Pmax    = MSRmax[Assetii] * 1.1
ntime   = length(profile)
segments= 52
factor  = floor(ntime/segments)
nvars   = factor

##Calculate optimal battery sizes
print("--Setting up optimization--")
# Pload = c(6,0,5,2,3,1,5,5,5,5,3,2,1,1,1,1,1,1,1)
# Pmax = 3
obj   = c(rep(0,nvars),1) #Objective function 

#Set up constraints (In a for loop, exclude this part, it has only to be done once)
sign = c(rep("<=",nvars),rep(">=",nvars),rep("<=",nvars),rep(">=",nvars))

Tmat = matrix(0,nvars,nvars)
Tmat[lower.tri(Tmat,diag=TRUE)] = 1
Amat = rbind(Tmat,Tmat,diag(-rep(1,nvars)),diag(-rep(1,nvars)))
Emat = c(rep(-1,nvars),rep(0,3*nvars)) #Set up storage variables
Amat = cbind(Amat,Emat)
bounds = list(lower = list(ind = c(1:(nvars+1)), val = c(rep(-Inf,nvars+1))),
              upper = list(ind = c(1:(nvars+1)), val = c(rep(Inf,nvars+1))))

print("--Optimizing...--")
for(ii in c(612,17,391)){
# for(ii in c(391)){ 
cl<-makeCluster(8) 
registerDoSNOW(cl)
tic()
Bsize = foreach(kk = 1:nMSR, .combine=rbind, .packages=c('Rglpk'), .verbose=FALSE) %dopar% {
   Assetii = kk
   profile = baseloadperMSR[Assetii,] + AllKVtechprofiles %*% KVScenariosperMSR[Assetii,,ii] + AllGVtechprofiles %*% GVScenariosperMSR[Assetii,,ii] 
   profile = profile[-((length(profile)-96-1):length(profile))]
   Pmax    = MSRmax[Assetii] * 1.1
   
   B = 0
   E = 0
   loadcycles = 0
   if((max(abs(profile))>Pmax & Pmax >0)){
      #Calculate battery size
      for(jj in 1:segments){
         Pload = profile[((jj-1)*factor+1):(jj*factor)]
         #Setup linear programming problem
         bmat = matrix(c(rep(0,2*nvars),(rep(Pmax,nvars)-Pload),(rep(-Pmax,nvars)-Pload)),4*nvars,1)
         
         #Run optimization
         
         
         A = Rglpk_solve_LP(obj, Amat, sign, bmat, bounds, max = FALSE)
         
         
         B = max(B,A$optimum)
         E = sum(A$solution[A$solution>0]) + E #Total stored energy ~ load cycles
          
#          #Calculate number of load cycles
#          nPeak = length(findValleys(A$solution, thresh=0)) + nPeak
      }
      loadcycles = E/B
   }
   return(c(B,E,loadcycles))
}
stopCluster(cl)
toc()

#Save results
print("--Saving results--")
Bsize[Bsize[,1]>0,] #Which MSR do have batteries?
table(Bsize[,1]>0&Bsize[,1]<100) #In how many cases is storage feasible (<100kWh)
plot(sort(Bsize[Bsize[,1]>0,1])[1:400])

save("Bsize",file='Batterysizes.Rdata')

### Generate CSV output
dt = data.table(1:nMSR,    # MSR number
                as.character(Bsize[,1]), # Necessary MSR battery size (kWh)
                as.character(Bsize[,2]), # Energy entered the battery in a year (kW)
                as.character(Bsize[,3])) # A/B ,number of load cycles
dt[dt==Inf] = NA
setnames(dt,c("MSR number",
              "Necessary MSR battery size (kWh)",
              "Energy entered the battery in a year (kW)",
              "A/B ,number of load cycles"))
filename = paste("01 MSR battery 2030 scenarionr ",ii,".csv")
write.csv2(dt,filename,row.names=FALSE)
}

# #Export an example battery profile
# C = A$solution[1:(length(A$solution)-1)]
# D = Pload[1:(length(A$solution)-1)]
# E = D-C
# dt = data.table(as.character(C), 
#                 as.character(D), 
#                 as.character(E), 
#                 as.character(1:(length(A$solution)-1)))
# setnames(dt,c("Batterij load (kW)",
#               "Vraag (kW)",
#               "MSR load door batterij (A-B) (kW)",
#               "Kwartiernummer"))
# filename = paste("08 MSR battery profile example.csv")
# write.csv2(dt,filename,row.names=FALSE)

#Get number of EANs per MSR
A = matprod_simple_triplet_matrix(KVEANtoMSR,matrix(rep(1,562370),562370,1))
dt = data.table(as.character(A))
setnames(dt,c("Number of EANs per MSR"))
filename = paste("10 EANs per MSR.csv")
write.csv2(dt,filename,row.names=FALSE)


# #Setup linear programming problem (small-scale example)
# Pload = c(6,0,5)
# MSRmax = 3
# nvars = length(Pload)
# obj = c(rep(0,nvars),1)
# 
# #Ax<=b
# sign = c(rep("<=",nvars),rep(">=",nvars),rep("<=",nvars),rep(">=",nvars))
# 
# 
# Amat = matrix(c(1 ,0 ,0 ,-1,
#                 1 ,1 ,0 ,-1,
#                 1 ,1 ,1 ,-1,
#                 1 ,0 ,0 ,0,
#                 1 ,1 ,0 ,0,
#                 1 ,1 ,1 ,0,
#                 -1 ,0 ,0 ,0,
#                 0 ,-1 ,0 ,0,
#                 0 ,0 ,-1 ,0,
#                 -1 ,0 ,0 ,0,
#                 0 , -1,0, 0,
#                 0 , 0, -1,0)
#               ,nrow =4*nvars,byrow = TRUE)
# 
# bmat = matrix(c(rep(0,2*nvars),(rep(MSRmax,nvars)-Pload),(rep(-MSRmax,nvars)-Pload)),4*nvars,1)
# bounds = list(lower = list(ind = c(1:(nvars+1)), val = c(rep(-Inf,nvars+1))),
#               upper = list(ind = c(1:(nvars+1)), val = c(rep(Inf,nvars+1))))
# 
# Rglpk_solve_LP(obj, Amat, sign, bmat, bounds, max = FALSE)

# i=c()
# j=c()
# v=c()
# progressbar = txtProgressBar(min = 0, max = nvars, initial = 0, char = "=", style = 3)
# for(ii in 1:nvars){
#   A = ii:nvars
#   setTxtProgressBar(progressbar,ii)
#   i             = c(i,A)
#   j             = c(j,rep(ii,nvars-ii+1))
#   v             = c(v,rep(1,nvars-ii+1))
# }
# close(progressbar)
# 
# Tmat = simple_triplet_matrix(i, j, v, nvars, nvars, dimnames = NULL) 