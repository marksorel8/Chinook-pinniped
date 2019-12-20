library(here)
library(TMB)
library(mvtnorm)
library(foreach)
library(doParallel)
 #------------------------------------------------------
# Functions
#------------------------------------------------------

#Function to process and load data as well as loading some helper functions
load_dat<-function(){
  #load helper functions
source(here("src","Chinook_pinn_mod_funcs_final.R")) 
out<-dataProcess(start_day = 70,end_day = 200,end_bon=209)    
save(out,file=here("pinn_out.Rdata"))
return(out)
}

#--------------------------------------------------
#--------------------------------------------------
#Function to create bounds for parameters
#--------------------------------------------------

make_bounds<-function(out){
#doesn't actually bound the random effects.
upper=c(rep(3,length(out$parameters$log_sigma)), #log st. dev. for penalized-complexity priors
        4,3, #logit lambda and log alpha
        rep(3,length(out$parameters$B)), # coefficients in travel time model
        6, #phi intercept
        rep(3,length(out$parameters$phi_betas)), #coefficients in survival model
        rep(6,length(out$parameters$A_mu)),  # poplation log mean Astoria departure
        rep(6,length(out$parameters$A_sigma)), # populaiton log st. dev. Astoria depature
        rep(1,length(out$parameters$A_year_beta)))#,
       # a_hypers=c(4,4,4,4)) # year effects on Astoria departure 
lower=c(rep(-7,length(upper)))
return(list(lower=lower,upper=upper))
}


#--------------------------------------------------
#--------------------------------------------------
#Function to fit model
#--------------------------------------------------

fit_mod<-function(data,params,fit_pop_Ast,fit_HMM,map_pars,haz_map,surv_map,fit){
  setwd(here("src"))
compile("pop_surv_final.cpp")
dyn.load(dynlib("pop_surv_final"))
  
# "flags" for which components of the full model to fit
data$fit_pop_Ast<-fit_pop_Ast #fit the populations-specific Astoria timing
data$fit_HMM<-fit_HMM   # fit the HMM


#creat a "map" which fixes  paramatersThis is done by assinging them "NA" factor values in the map.
map<-create_Map(map_pars,params)

if(!is.null(haz_map)){
#Map hazard model of transitiona proability (travel time)
map$B <-factor(haz_map) 

params$B[which(is.na(map$B))]<-0 #fix mapped coefs at 0.0.
}

if(!is.null(surv_map)){
#Map coefficients in daily survival logit model
map$phi_betas<-factor(surv_map) #order is 1)sea lions, 2)adipose clip, 3) temperature

params$phi_betas[which(is.na(map$phi_betas))]<-0 ##fix mapped coefs at 0.0.
}

#sigmas for penalized complexity normal priors
map$log_sigma<-1:length(params$log_sigma)
map$log_sigma[is.na(c(map$B,map$phi_betas))]<-NA #Map out (fix) sigmas for paramaters that are mapped out
map$log_sigma<-factor(map$log_sigma)

start_time <- Sys.time()
#Initialize model
mod<-MakeADFun(data,params,
                    random = c("phi_betas","B"),
                    DLL="pop_surv_final",silent=T,map=map,
                    inner.control = list(tol10=0))


#optimize model
if(fit){ mod_fit<-nlminb(mod$par, mod$fn, mod$gr,
                     control = list(iter.max=1e3,eval.max=1e3))}
end_time <- Sys.time()


print(end_time - start_time)

return(mod)
}

#--------------------------------------------------
#--------------------------------------------------
#draw parameter sets from a Multivariate normal for the HMM
#--------------------------------------------------

draw_params<-function(mod){
sd_mod<-sdreport(mod)

sum_sd<-summary(sd_mod)[1:15,]   # standard deviations of paramaters
upper<-qnorm(.975,sum_sd[,1],sum_sd[,2])
lower<-qnorm(.025,sum_sd[,1],sum_sd[,2])
sum_sd<-cbind(sum_sd,lower,upper) # add CI
print(sum_sd) #view

set.seed(1)
par_draws<-rmvnorm(1000,mod$env$last.par.best,sd_mod$cov[-which(sd_mod$cov[1,]==0),
                                                               -which(sd_mod$cov[1,]==0)])

return(par_draws)
}


#--------------------------------------------------
#--------------------------------------------------
#Function to fit the population-specific timing model for each parameter draw
#--------------------------------------------------

fit_pop_time<-function(mod,data,par_draws,Niter=1000){

ncores<-detectCores()-2    # number of cores to use
  
cl <- makeCluster(ncores)  #settup clusters for parallel computing
registerDoParallel(cl)

pars<-fixedParams(mod,1:6) #declare starting paramater list to fill with each parameter set within loop

mod_map<-mod$env$map

#loop through parameter sets and fit Astoria departure (migration) timing and store fits
#this takes about three hours
system.time(
  results<-foreach(i =iter(1:Niter,chunksize=Niter/ncores),.packages=c('TMB','mvtnorm','here'),.export = c("fit_mod","create_Map"),.multicombine=TRUE,.inorder=FALSE) %dopar%{

    pars$log_sigma[!is.na( mod_map$log_sigma)]<-par_draws[i,colnames(par_draws)=="log_sigma"]
    pars$logit_lambda<-par_draws[i,colnames(par_draws)=="logit_lambda"]
    pars$log_alpha<-par_draws[i,colnames(par_draws)=="log_alpha"]
    pars$B[!is.na( mod_map$B)]<-par_draws[i,colnames(par_draws)=="B"]
    pars$phi_intercept<-par_draws[i,colnames(par_draws)=="phi_intercept"]
    pars$phi_betas[!is.na( mod_map$phi_betas)]<-par_draws[i,colnames(par_draws)=="phi_betas"]
  

fit1<-fit_mod(data,pars,fit_pop_Ast=1,fit_HMM=0,map_pars=1:6,haz_map=NULL,surv_map=NULL,fit=TRUE)


   #calculate covariance matrix
    sdRep<- sdreport(fit1)
    #draw parameter sets from multivariate normal based on MLE
    set.seed(1)
    rmvnorm(10,fit1$env$last.par.best,sdRep$cov.fixed)
  }#end of loop
  
)



stopCluster(cl) #shut down parallel computing clusters
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister()


results2<-array(unlist(results),dim=c(dim(results[[1]]),Niter))
save(results2,file=here("AvePin_3par.RData")) #save parameter sets from population-specific departure model
return(results2)
}

#--------------------------------------------------
#--------------------------------------------------
#loop through the paramater sets and calculate population specific survival and the mortality attributable to increased pinniped abundance.
#--------------------------------------------------

pop_surv<-function(par_draws,pop_time_boot,Niter,...){
results3<-matrix(NA,nrow=180*Niter,ncol=7)#matrix to hold results

full_mod<-fit_mod(data=out$data,params=out$parameters,fit_pop_Ast=1,fit_HMM=1,map_pars=NULL,haz_map=NULL,surv_map=NULL,fit=FALSE)

#this takes about 20 minutes
for ( i in 1:Niter){#loop through paramater sets from the HMM
  for (j in 1:10){#loop through paramater sets for the departure timing model for each HMM paramter set
    
    mod_iter_rep<-full_mod$report(c(par_draws[i,],pop_time_boot[j,,i]))#report with the paramater sets
    results3[((i*180 -179)+((j-1)*18)):((i*180 -179)+((j)*18)-1),1:6]<-mod_iter_rep$pop_yr_mort[,1:6]#save population specific mortality
    results3[((i*180 -179)+((j-1)*18)):((i*180 -179)+((j)*18)-1),7]<-(1-mod_iter_rep$pop_yr_mort[,8])/(1-mod_iter_rep$pop_yr_mort[,7])#calculate and store population-specific percent change in survival
    if(((i*j)/Niter)%%1==0) print(i)#keep track of progress
  }
}

save(results3,file=here("boot_pop_surv.Rdata"))# save output

return(results3)
}


##
summarize_change_func<-function(pop_surv_out){
  #matrix of population specific percent changes in survival accross year-groups for each parameter set
  per_change<-matrix(pop_surv_out[,7],ncol=18,byrow = TRUE)
  
  print(c("range of population-specific changes",range(apply((per_change), 2, median))))
  
  #quantiles of percent change across parameter sets for each population
  per_change_qant<-apply(per_change, 2, quantile,probs=c(.025,.25,.5,.75,.975))
  
  
  colnames(per_change_qant)<-levels(out$dat[[2]]$Pop)
  
  #mortality
  print(1-per_change_qant)
  
  #look at the distribution of median population specific mortalities to split them into the highest and lowest (earliest and latest)
  hist(1-per_change_qant[3,],breaks=18)
  
  #percent change grouped by early and late
  Late<-1-per_change_qant[,(1-per_change_qant[3,])<=.15]
  Early<-1-per_change_qant[,(1-per_change_qant[3,])>=.15]
  
  print("Average early pops")
  print(rowMeans(Early))
  print("average late pops")
  print(rowMeans(Late))
  
}




#--------------------------------------------------
#  End functions
#--------------------------------------------------
#  Start analysis
#-------------------------------------------------


#---------------------------
#load data and paramater inits
set.seed(1234)
 out<-load_dat()
#---------------------------
#bounds on parameters
bounds<-make_bounds(out)
#---------------------------
#Function to fit model
HMM_only<-fit_mod(data=out$data,params=out$parameters,fit_pop_Ast=0,fit_HMM=1,map_pars=7:9,haz_map=1:3,surv_map=1:3,fit=TRUE)

 HMM_only$env$last.par.best
 HMM_only$par
 HMM_only$fn()
 HMM_only$gr()
 assessFunc(HMM_only, out$dat,70,200,209)
 
 
 
 new_pin_dat_mod<-fit_mod(data=out$data,params=fixedParams(HMM_only,1:6),fit_pop_Ast=1,fit_HMM=0,map_pars=1:6,haz_map=NULL,surv_map=NULL,fit=FALSE)
 
 
 sim1<-new_pin_dat_mod$simulate()
 sim2<-sim1$sim
 hist(sim2)
sim_surv<-data.frame(surv=as.numeric(sim2!=140),TT_bon=ifelse(as.numeric(sim2!=140),sim2-out$data$relDOY,NA))
write.csv(sim_surv,file=here("data","sim_surv.csv"))

 
#---------------------------
#draw parameter sets from a Multivariate normal for the HMM
par_draws<-draw_params(HMM_only)

 start1<-Sys.time()
#fit the population-specific timing model for each HMM parameter draw
if(file.exists(here("AvePin_3par.RData"))){load(here("AvePin_3par.RData"))}else{
pop_time_boot<-fit_pop_time(mod=HMM_only,data=out$data,par_draws=par_draws ,Niter=1000)
save(pop_time_boot,file=here("AvePin_3par.RData"))

}
 end1<-Sys.time()
#---------------------------
#loop through the full paramater sets and calculate population specific survival and the mortality attributable to increased pinniped abundance.
 start2<-Sys.time()
if(file.exists(here("boot_pop_surv.Rdata"))){load(here("boot_pop_surv.Rdata"))}else{
pop_surv_out<-pop_surv(par_draws=par_draws ,pop_time_boot=pop_time_boot,Niter = 1000)
save(pop_surv_out,file=here("boot_pop_surv.Rdata"))
}
#---------------------------
 end2<-Sys.time()
#plotting
png(file = here("pop_surv.png"),units="in",res=300, width=4,height=4)
plot_pop_surv(pop_surv_out,start_day = 70, end_day = 200)
dev.off()
#---------------------------

#Fit population timing model with the paramaters in the HMM set at the MLE
pop_time_fit<-fit_mod(data=out$data,params=fixedParams(HMM_only,1:6),fit_pop_Ast=1,fit_HMM=0,map_pars=1:6,haz_map=NULL,surv_map=NULL,fit=TRUE)

assessFunc(pop_time_fit, out$dat,70,200,209)

png(file = here("change_plot.png"),units="in",res=300, width=6,height=5)
coef_plot(report= pop_time_fit$env$last.par.best,
          start_day=70,surv_change=pop_surv_out)
dev.off()

#---------------------------

summarize_change_func(pop_surv_out)

# end of script



#simulation
sim<-matrix(NA,50,length(HMM_only$simulate()$sim))
for ( i in 1:50){
  sim[i,]<-HMM_only$simulate()$sim
}
hist(sim)

sim_result<-matrix(NA,50,length(HMM_only$env$last.par.best))
for ( i in 1:50){
  dat<-out$data
  dat$bonDOY_unk_pop<-sim[i,]
  HMM_only2<-fit_mod(data=dat,params=out$parameters,fit_pop_Ast=0,fit_HMM=1,map_pars=7:9,haz_map=1:3,surv_map=1:3,fit=TRUE)
  sim_result[i,]<- HMM_only2$env$last.par.best
} 



for ( i in 1:ncol(sim_result)){
  hist(sim_result[,i],main=names(HMM_only$env$last.par.best)[i])
  abline(v=HMM_only$env$last.par.best[i],col="red")
}
apply(sim_result, 2, mean)
sim_fits