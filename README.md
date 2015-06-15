Analytically solving a soil C cycling model using the Newton-Raphson method. Implemented using the `stode` function from the `rootSolve` package in R.
-----
This Repository is an example of how to use the `stode` funtion in the `rootSolve` package in R, to solve for model outputs at steady-state. I use the "carbon only" model from Schimel and Weintraub 2003 as the example model.


Model Code
-----
This is how I code the model and normally run to steady state. 

```{r}
t= 800 # number of days to run the model for
day = 0 #start at day 0

#initial pools
MBC = 28 #gleaned from figures in the manuscript
EnzC= 1.4 #multiplied initial MBC by Ke, 0.05.

#parameters
Kd = 1 #this is the K and SOC values collapsed as 1 term
Ke = 0.05 #this is the fraction of Uc that goes to enzymes
Kl = 0.05 #this is the enzyme turnover rate per day
Km = 0.022 #microbial maintenance rate for biomass
SUE = 0.5 #taken from manuscript- they say this is from literature
Kes = 0.3 # half saturation of enzymes on substrate constance
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
  #end model loop
}
out<-data.frame(out)
```
testing figures
```{r}
hist(iris[[2]])
```


Once the model has run 800 days the pools are near steady state.
```{r, echo=FALSE}
par(mfrow=c(1,2))
plot(MBC~day,data=out,xlab='day',ylab='microbial biomass (mg / g soil)')
plot(EnzC~day,data=out,xlab='day',ylab='enzyme carbon (mg / g soil)')
```


