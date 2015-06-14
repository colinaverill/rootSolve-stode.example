#testing rootSolve with the schimel and weintraub 2003 C-only model. 

#initial pools
MBC = 28 #gleaned from figures in the manuscript
EnzC= 1.4 #multiplied initial MBC by Ke, 0.05.
day = 0 #start at day 0

#parameters
Kd = 1 #this is the K and SOC values collapsed as 1 term
Ke = 0.05 #this is the fraction of Uc that goes to enzymes
Kl = 0.05 #this is the enzyme turnover rate per day
Km = 0.022 #microbial maintenance rate for biomass
SUE = 0.5 #taken from manuscript- they say this is from literature
Kes = 0.3 # half saturation of enzymes on substrate constance
t= 800 # number of days to run the model for
params<-c(Kd,Ke,Kl,Km,SUE,Kes,t)

#choose model outputs you want to track.need to specify again at end of model code (how do I change this so you only have to specify it here?)
outputs<-c('day','MBC','EnzC','R')
#create a matrix to save model outputs 
out<- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))


#run simulation t times
for(i in 1:t){
  #fluxes!
  Dc = Kd * (EnzC / (Kes + EnzC)) #decomposition
  Uc = Dc #upatke of C
  EPc = Ke * Uc #enzyme production
  ELc = Kl * EnzC # enzyme loss 
  Re = EPc * ((1-SUE)/SUE) #enzyme production respiration
  Rm = Km * MBC #maintenance respiration
  Rg = (Uc - EPc - Re - Rm)*(1 - SUE) #growth respiration
  MBCg = (Uc - EPc - Re - Rm - Rg) #biomass growth
  R = Re + Rm + Rg
  #pools!
  MBC = MBC + MBCg 
  EnzC = EnzC + EPc - ELc
  day = day + 1
  #update output matrix
  out[i,]<- c(day,MBC,EnzC,R)
}
out<-data.frame(out)

require(rootSolve)

#model as a function
#the model no longer updates pool sizes and iterates through time. It is just a function. At the end we specify the differential equations (dMBC and dEnz), the model will try to solve for (it iterates until the values of these functions is 0)
model<-function(t,y,pars){
  with(as.list(c(y,pars)),{
    #fluxes!
    Dc = Kd * (EnzC / (Kes + EnzC)) #decomposition
    EPc = Ke * Dc #enzyme production
    ELc = Kl * EnzC # enzyme loss 
    Re = EPc * ((1-SUE)/SUE) #enzyme production respiration
    Rm = Km * MBC #maintenance respiration
    Rg = (Uc - EPc - Re - Rm)*(1 - SUE) #growth respiration
    MBCg = (Uc - EPc - Re - Rm - Rg) #biomass growth
    R = Re + Rm + Rg
    #pools to solve for!
    dMBC = MBCg
    dEnzC = EPc - ELc
    list(c(dMBC,dEnzC))
  })
}

#specify parameter vector for stode
pars<-c(Kd=Kd,Ke=Ke,Kl=Kl,Km=Km,SUE=SUE,Kes=Kes)

#specify initial pool sizes, you are solving for steady state of these values.
y<-c(MBC=28,EnzC=1.4)

ST<- stode(y=y,func=model,parms=pars,pos=T)