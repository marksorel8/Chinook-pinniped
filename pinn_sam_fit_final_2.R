library(here)
library(TMB)
library(TMBhelper)
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
mod_fit<-NULL
if(fit){ mod_fit<-TMBhelper::fit_tmb(mod, mod$fn, mod$gr,newtonsteps = 1,getReportCovariance=TRUE)}
  
end_time <- Sys.time()


print(end_time - start_time)

return(list(mod=mod,mod_fit=mod_fit))
}

#--------------------------------------------------
#--------------------------------------------------
#Function to fit the population-specific timing model for each parameter draw
#--------------------------------------------------

fit_pop_time<-function(mod,data,par_draws,Niter=1000,ncores=1){

cl <- makeCluster(ncores)  #settup clusters for parallel computing
registerDoParallel(cl)

pars<-fixedParams(mod,1:6) #declare starting paramater list to fill with each parameter set within loop

mod_map<-mod$env$map

#loop through parameter sets and fit Astoria departure (migration) timing and store fits
#this takes about three hours

  results<-foreach(i =iter(1:Niter,chunksize=Niter/ncores),.packages=c('TMB','mvtnorm','here'),.export = c("fit_mod","create_Map"),.multicombine=TRUE,.inorder=FALSE) %dopar%{


    pars$log_lambda<-par_draws[i,colnames(par_draws)=="log_lambda"]
    pars$log_alpha<-par_draws[i,colnames(par_draws)=="log_alpha"]
    pars$B[!is.na( mod_map$B)]<-par_draws[i,colnames(par_draws)=="B"]
    pars$phi_intercept<-par_draws[i,colnames(par_draws)=="phi_intercept"]
    pars$phi_betas[!is.na( mod_map$phi_betas)]<-par_draws[i,colnames(par_draws)=="phi_betas"]
  

fit1<-fit_mod(data,pars,fit_pop_Ast=1,fit_HMM=0,map_pars=1:6,haz_map=NULL,surv_map=NULL,fit=TRUE)



    #draw parameter sets from multivariate normal based on MLE
    set.seed(1234)
    rmvnorm(10,fit1$mod_fit$par,fit1$mod_fit$SD$cov.fixed)
  }#end of loop
  




stopCluster(cl) #shut down parallel computing clusters
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister()


results2<-array(unlist(results),dim=c(dim(results[[1]]),Niter))
return(results2)
}

#--------------------------------------------------
#--------------------------------------------------
#loop through the paramater sets and calculate population specific survival and the mortality attributable to increased pinniped abundance.
#--------------------------------------------------

pop_surv<-function(par_draws,pop_time_boot,Niter,...){


results3<-array(NA,dim=c(18,8,Niter*10))
  
full_mod<-fit_mod(data=out$data,params=out$parameters,fit_pop_Ast=1,fit_HMM=1,map_pars=NULL,haz_map=NULL,surv_map=NULL,fit=FALSE)

#this takes about 20 minutes

for ( i in 1:Niter){#loop through paramater sets from the HMM
  for (j in 1:10){#loop through paramater sets for the departure timing model for each HMM paramter set
  
    mod_iter_rep<-full_mod$mod$report(c(rep(0,6),par_draws[i,],pop_time_boot[j,,i]))#report with the paramater sets
    results3[,,((i-1)*10)+j]<-mod_iter_rep$pop_yr_mort

  }
}

save(results3,file=here("boot_pop_surv.Rdata"))# save output

return(results3)
}


##
summarize_change_func<-function(pop_surv_out){
  #matrix of population specific percent changes in survival accross year-groups for each parameter set
  
  surv<-function(column,time_vec=NULL){
     x<-apply(1-pop_surv_out[,column,],1,quantile,probs=c(.025,.25,.5,.75,.975))
  colnames(x)<-levels(out$dat[[2]]$Pop)
  if(is.null(time_vec))time_vec<-(x[3,])<=.83
  
  Late<-x[,time_vec]
  Early<-x[,!time_vec]
  print("Average late pops")
  print(rowMeans(Early))
  print(colnames(Early))
  print("average early pops")
  print(rowMeans(Late))
  print(colnames(Late))
  return (time_vec)
  }
  
  print("baseline survival")
  baseline_surv<-surv(7)
  
  print("New survival")
  new_surv<-surv(8,time_vec=baseline_surv)

  
  
  mort<-(pop_surv_out[,8,]-pop_surv_out[,7,])
  
  print(c("range of population-specific changes",range(apply((mort), 2, median))))
  
  #quantiles of percent change across parameter sets for each population
  mort_qant<-apply(mort, 1, quantile,probs=c(.025,.25,.5,.75,.975))
  
  
  colnames(mort_qant)<-levels(out$dat[[2]]$Pop)
  

  
  #look at the distribution of median population specific mortalities to split them into the highest and lowest (earliest and latest)
  hist(mort_qant[3,],breaks=18)
  
  #percent change grouped by early and late
  Late<-mort_qant[,(mort_qant[3,])<=.15]
  Early<-mort_qant[,(mort_qant[3,])>=.15]
  print("additional mortality")
  print("Average early pops")
  print(rowMeans(Early))
  print(colnames(Early))
  print("average late pops")
  print(rowMeans(Late))
  print(colnames(Late))
  
  return(mort)
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

out$data$ad_mort <-1
 out$data$ad_mig <-1
#fit HMM model only
 HMM_only<-fit_mod(data=out$data,params=out$parameters,fit_pop_Ast=0,fit_HMM=1,map_pars=7:10,haz_map=1:3,surv_map=1:3,fit=TRUE)
 
 dev.new()
 assessFunc(HMM_only, out$dat,70,200,209, plot_TT = T,plot_surv = T,print_outs = FALSE,cex=.5,plot_out = F)


  #extra mortality for hatchery fish
 ast_tag_dat<-out$data #list of data with tag day, tag year, and origin (hatchery/wild)
 all_zeros_prob<-HMM_only$mod$report()$all_zeros_prob # mortality 
  harvest_mort_est<-numeric(sum(ast_tag_dat$relOrigin)) #empty vector for all hatchery fish tagged in study
 for ( i in 1:length(harvest_mort_est)){#loop over hatcheryfish in study
   n<-which(ast_tag_dat$relOrigin==1)[i]
   harvest_mort_est[i]<- all_zeros_prob[[2]][ast_tag_dat$relDOY[i]+1,ast_tag_dat$rel_Year[i]+1]-all_zeros_prob[[1]][ast_tag_dat$relDOY[i]+1,ast_tag_dat$rel_Year[i]+1]# difference between modeled mortality of hatchery and wild
 }
  mean(harvest_mort_est)#average additional mortality experienced by hatchery fish in study
  
 #---------------------------
#draw parameter sets from a Multivariate normal for the HMM
  set.seed(1234)
par_draws<- mvtnorm::rmvnorm(1000,HMM_only$mod_fit$SD$value[1:9],HMM_only$mod_fit$SD$cov[1:9,1:9])
colnames(par_draws)<-names(HMM_only$mod_fit$SD$value[1:9])
  

#fit the population-specific timing model for each HMM parameter draw
if(file.exists(here("results","AstPin_3par2_05242000.RData"))){load(here("results","AstPin_3par2_05242000.RData"))}else{

out$data$ad_mort<-0 #don't incldue mortality and travel time parameters in adreport
out$data$ad_mig <-0

start_time<-Sys.time()
pop_time_boot2<-fit_pop_time(mod=HMM_only$mod,data=out$data,par_draws=par_draws ,Niter=1000,ncores = (detectCores()-2))
save(pop_time_boot,file=here("AstPin_3par2.RData"))
Sys.time()-start_time
}

 #---------------------------
#loop through the full paramater sets and calculate population specific survival and the mortality attributable to increased pinniped abundance.
 start2<-Sys.time()
if(file.exists(here("results","boot_pop_surv_05242000.Rdata"))){load(here("results","boot_pop_surv_05242000.Rdata"))}else{
pop_surv_out<-pop_surv(par_draws=par_draws ,pop_time_boot=pop_time_boot2,Niter = 1000)
save(pop_surv_out,file=here("boot_pop_surv.Rdata"))
}
 end2<-Sys.time()
end2-start2

mort_pop_mat<-summarize_change_func(pop_surv_out)


#---------------------------

#plotting


#Fit population timing model with the paramaters in the HMM set at the MLE
pop_time_fit<-fit_mod(data=out$data,params=fixedParams(HMM_only$mod,1:6),fit_pop_Ast=1,fit_HMM=0,map_pars=1:6,haz_map=NULL,surv_map=NULL,fit=TRUE)

pop_time_rep<-pop_time_fit$mod$report()

dev.new()
assessFunc(pop_time_fit, out$dat,70,200,209,plot_pop_Ast = T, plot_TT = F,plot_surv = F,print_outs = FALSE,cex=.5,plot_out = F)


png(file = here("results","Astoria time.png"),units="in",res=300, width=4,height=3.5)
coef_plot(report= pop_time_fit$mod$env$last.par.best,
          start_day=70,mort=mort_pop_mat,plot_which=1)
dev.off()


png(file = here("results","change_plot.png"),units="in",res=300, width=4,height=3.5)
coef_plot(report= pop_time_fit$mod$env$last.par.best,
          start_day=70,mort=mort_pop_mat,plot_which=2)
dev.off()



####### plot astoria pinipeds
png("astoria_CLS.png", width=4,height=3.5,units="in",res=300)
#dev.new(width=4,height=4,units="in",res=300)
par(mfrow=c(3,2),mar=rep(0,4),oma=c(5,10,2,2),cex=.35)
for ( i in 1:6){
  
  plot(1,type="n",ylim=exp(range(out$pinn$log_pinn_Mat[1:139,1:6])),xlim=c(1,139),main="",xlab="",ylab="",yaxt="n",xaxt="n")
  if(((i/2)%%1)>0) axis(2,cex.axis=2.7)
  points(1:139,exp(out$pinn$log_pinn_Mat[1:139,i]),type ="l",lwd=4)
  if(i>4)axis(1,at=c(91,121,152,182)-70,labels = c("1-Apr","","1-Jun",""),cex.axis=2.7,padj=1)
  text(x=(139/2)+10,y=exp(max(out$pinn$log_pinn_Mat[1:139,1:6]))*0.75,labels = i+2009, pos=3,cex=2.7,font=2)
}
mtext("# Sea lions",2,5,T,cex=1.25)
dev.off()


# end of script
