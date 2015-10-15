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
library(Rglpk) 
library(data.table)
library(ggplot2)

#Set number of timesteps and profile
nvars = 10
Pload = c(15,0,0,0,5,6,7,8,9,0)
Pmax    = 5

#Setup Objective function 
obj   = c(rep(0,nvars),1) 

##Setup Constraints
#Setup signs
sign = c(rep("<=",nvars),rep(">=",nvars),rep("<=",nvars),rep(">=",nvars))

#Setup linear inequality constraints
Tmat = matrix(0,nvars,nvars)
Tmat[lower.tri(Tmat,diag=TRUE)] = 1
Amat = rbind(Tmat,Tmat,diag(-rep(1,nvars)),diag(-rep(1,nvars)))
Emat = c(rep(-1,nvars),rep(0,3*nvars))    #Set up storage variables
Amat = cbind(Amat,Emat)
bmat = matrix(c(rep(0,2*nvars),(rep(Pmax,nvars)-Pload),(rep(-Pmax,nvars)-Pload)),4*nvars,1)

#Setup bounds
bounds = list(lower = list(ind = c(1:(nvars+1)), val = c(rep(-Inf,nvars+1))),
              upper = list(ind = c(1:(nvars+1)), val = c(rep(Inf,nvars+1))))

#Run optimization
A = Rglpk_solve_LP() #Add the correct arguments

#Plot resuls
X = data.frame(time = 1:nvars,P_use=A$solution[1:nvars])
ggplot(X,aes(x=time,y=P_use))+geom_line()
print(paste('The optimal battery size is:',A$solution[nvars+1],'kWh'))
