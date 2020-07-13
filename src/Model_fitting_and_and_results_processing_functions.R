# a colloection of functions for conducting the analysis of population-specific mortlaity associated with increased sea lion abundance in 2013-2105 relative to 2010-2012

#Mark Sorel 7-12-2020

#--------------------------------------------------
#--------------------------------------------------
#Function to process and load data as well as loading some helper functions
#--------------------------------------------------
#returns a list containing a list of model data and a list of initial parameters
load_dat<-function(){
  #load helper functions
  source(here("src","data_processing.R")) 
  out<-dataProcess(start_day = 70,end_day = 200,end_bon=209)    
  save(out,file=here("pinn_out.Rdata"))
  return(out)
}

#--------------------------------------------------
#--------------------------------------------------
#Function to create bounds for parameters
#--------------------------------------------------
#takes a list with data and paramaters from the "load_dat" function. Returns a list of uper and lower bounds for parameters
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
#helper function to create TMB map that "maps out" entire paramater vectors
#--------------------------------------------------
#takes vector of indices of paramaters to be mapped out, and the list of paramater initial values. returns a list that is a "map" of which paramaters not to optimize. 
create_Map<-function(par_nums,params){
  
  map<-list()
  for ( i in par_nums){
    map[[names(params[i])]]<-factor(rep(NA,length(params[[i]])))
  }
  return(map)
}


#--------------------------------------------------
#--------------------------------------------------
#Function to fit model
#--------------------------------------------------
#takes lists of data and initial paramater values, and flags for whether to optimize the populaiton specific Astoria departure time sub model (fit_pop) and/or the hiddem markvov model (fit_HMM), as well as vectors of parameters to "map out" (i.e., not fit) in the survival model(map_pars) and the daily dam passage probability model (Haz_map), and whether to conduct any optimization at all (fit). Returns a fitted model (mod) and a summary or output (mod_fit)
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
#function to format parameter values for input based on a model fit. 
#--------------------------------------------------
#takes a fitted model object and vector of parameter indices to fill with values from the fitted model object. Returns a list of parmater values with the selected values filled in from the fitted model object
fixedParams<-function(mod,pars=1:6){
  parameters2<-out$parameters
  parNames=names(parameters2)[pars]
  for(i in parNames){
    
    parameters2[[i]][which(!is.na( mod$env$map[[i]]))]<-mod$env$last.par.best[names(mod$env$last.par.best)==i]
    
    if(is.null( mod$env$map[[i]])){ parameters2[[i]]<-mod$env$last.par.best[names(mod$env$last.par.best)==i]}
    
    
  }
  parameters2
}


#--------------------------------------------------
#--------------------------------------------------
#Function to fit the population-specific timing model for each parameter draw
#--------------------------------------------------
#takes model (mod), data, paramaters for the HMM drawn from a multivariate normal distirbution representing the covariance matrix (par_draws), number of iterations, and number of cores to use. Returns an array of parameter values for the population- an d year specific Astoria Departure timing model, where 10 values are returned for each HMM paramater set (iteration), and the 10 values are drawn for a multivariate normal representing the covariance matrix of the population- and year-specific migration timing model fit with a given HMM paramater set. 
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
#loop through the paramater sets and calculate population specific survival 
#--------------------------------------------------
#takes arrays of paramater sets for the HMM (par_draws) and the population-specific Astoria departure models(pop_time_boot), as well as the number of paramater sets (Niter). Returns an array of population- and year-specific survival estimates.
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


#--------------------------------------------------
#loop through the paramater sets and calculate the mortality attributable to increased pinniped abundance.
#--------------------------------------------------
#takes an array of population- and year-specific survival estimates. Prints summary statistics for survival and mortality of earliest and latest migrating populations, and returns a matrix of estimates population-specific mortality associated with the increase in sea lions in 2013-2015 relative to 2010-2012, with an estimate for each paramater set.  
summarize_change_func<-function(pop_surv_out,plot=FALSE){
  #matrix of population specific percent changes in survival accross year-groups for each parameter set
  
  surv<-function(column,time_vec=NULL){
    x<-apply(1-pop_surv_out[,column,],1,quantile,probs=c(.025,.25,.5,.75,.975))
    colnames(x)<-levels(out$dat[[2]]$Pop)
    if(is.null(time_vec))time_vec<-(x[3,])<=.83
    
    Late<-x[,time_vec]
    Early<-x[,!time_vec]
    print(" late pops")
    print(rowMeans(Early))
    print(colnames(Early))
    print("")
    print(" early pops")
    print(rowMeans(Late))
    print(colnames(Late))
    print("")
    return (time_vec)
  }
  
  print("2010-2012 average survival")
  baseline_surv<-surv(7)
  
  print("2013-2015 average survival")
  new_surv<-surv(8,time_vec=baseline_surv)
  
  
  
  mort<-(pop_surv_out[,8,]-pop_surv_out[,7,]) #calculate mortality
  
 
  #quantiles of percent change across parameter sets for each population
  mort_qant<-apply(mort, 1, quantile,probs=c(.025,.25,.5,.75,.975))
  
  
  colnames(mort_qant)<-levels(out$dat[[2]]$Pop) 
  
  
  
  #look at the distribution of median population specific mortalities to split them into the highest and lowest (earliest and latest)
  if(plot) hist(mort_qant[3,],breaks=18)
  
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
