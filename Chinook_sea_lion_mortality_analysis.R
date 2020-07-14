# This script conduucts the analysis of the association between increased sea lion abundance in 2013-2015 and mortality of adult CHinook salmon in the Lower COlumbia River.

# Mark Sorel 7-12-20

library(here)
library(TMB)
library(TMBhelper)
library(mvtnorm)
library(foreach)
library(doParallel)

#load function for model fitting and results processing
source(here("src","Model_fitting_and_and_results_processing_functions.R"))

#load function for plotting results
source(here("src","plotting_functions.R"))

#---------------------------
#load data and paramater inits. Also loads functions for plotting
out<-load_dat()
#---------------------------
#bounds on parameters
bounds<-make_bounds(out)
#---------------------------
#beggin model fitting
out$data$ad_mort <-1  #set flad to calculate delta method standard errors for survival model
 out$data$ad_mig <-1  #set flag to calculate delta method standard errors for Bonneville passage probabilty (travel time) model
#fit HMM model only
 HMM_only<-fit_mod(data=out$data,params=out$parameters,fit_pop_Ast=0,fit_HMM=1,map_pars=7:10,haz_map=1:3,surv_map=1:3,fit=TRUE)
 
 assessFunc(HMM_only, out$dat,70,200,209, plot_TT = T,plot_surv = T,print_outs = TRUE,cex=.5,plot_out = F)


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
  

#fit the population-specific timing model for each HMM parameter draw (this takes hours to run. With 50 cores it only took 45 mintues however. The Rdata file is saved in the GitHub repository however.)
if(file.exists(here("results","AstPin_3par2_05242000.RData"))){load(here("results","AstPin_3par2_05242000.RData"))}else{

out$data$ad_mort<-0 #don't incldue mortality and travel time parameters in adreport
out$data$ad_mig <-0

start_time<-Sys.time()
pop_time_boot2<-fit_pop_time(mod=HMM_only$mod,data=out$data,par_draws=par_draws ,Niter=1000,ncores = (detectCores()-2))
save(pop_time_boot,file=here("AstPin_3par2.RData"))
Sys.time()-start_time
}

 #---------------------------
#loop through the full paramater sets and calculate population specific survival and the mortality attributable to increased pinniped abundance. (This takes about an hours to run)
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
