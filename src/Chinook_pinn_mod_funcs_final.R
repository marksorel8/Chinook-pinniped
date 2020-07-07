
###################################################
#              Data loading and wrangling
####################################################


dataProcess<-function(start_day=75,end_day=200,end_bon=210,ad_mort=1,ad_mig=0){ #function to read and process data, returns list of data for model, and a list of initial parameters for model, as well as data for plotting
  
  library(here)
  
  #read data with Astoria tagging, survival, and Bonneville detection date data
  dat<-read.csv(here("data","surv_dat_MS.csv"))
  
  #load Astoria sea lion count data
  pinnDat<-read.csv(here("data","pinnCounts.csv"))
  
  #create a time series with NAs for days with no pinniped counts 
  ##create sequence of zeros
  dateSeq<-paste(rep(2010:2015,each=365),rep(1:365,times=length(2010:2015)))
  ##make data frame
  pinnTS<-data.frame(year=rep(2010:2015,each=365),doy=rep(1:365,times=length(2010:2015)),
                     count=pinnDat$Sum.of.Max.of.Count[match(dateSeq,paste(pinnDat$year,pinnDat$julian))])
  ##make time series of sea lion counts
  pinnTS2<-ts(pinnTS$count,start=c(2010,1),frequency=365)
  
  #plot(pinnTS2)
  #fit a loess model to the pinniped counts to interpolate missing values.
  #This model could be integrated with the rest of the model I suppose 
  loessTest<-loess(log(pinnTS$count+1)~c(1:length(pinnTS$count)),span=.02)
  #make predictions for daily pinniped counts from loess model
  predPinn<-ts(predict((loessTest),
                       newdata =data.frame( x=c(1:length(pinnTS$count)))),
               start=c(2010,1),frequency = 365)
  
    #points(exp(predPinn),type="l",col="blue")

  #make a matrix of Astoria sea lion counts (with a column for each year) 
  log_pinn_Mat<-matrix((predPinn),ncol=6)[start_day:282,]
  
  
  #load bonneville pinniped counts [start_day:282,]
  bon_pinn_Mat<-source(here("src","bon_pin_ts_func.R"))$value
  
  #make it into one long time secies                     
  bon_pinn_Mat<-c(bon_pinn_Mat)
 
  #fit loess smoother to log +1
  log_bon_pinn_loess<-loess(log(bon_pinn_Mat+1)~c(1:length(bon_pinn_Mat)),span=.02)
  
  
  bonn_pinn_pred<-(predict((log_bon_pinn_loess),
                              newdata =data.frame( x=c(1:length(bon_pinn_Mat)))))
  

#plot(exp(bonn_pinn_pred),type="l")
#abline(v=rep(c(70,210),times=6)+rep((0:5)*365,each=2)) 
   #make a matrix of Bonneville sea lion counts (with a column for each year) 
  log_bon_pinn_Mat<-matrix((bonn_pinn_pred),ncol=6)[start_day:282,]
  
  
  #scale
  ##Because a covariate value will be used for each day of each fish's capture history, we need to scale the covariates based on how they are used in the model. This function creates a vector of a covariate equivalent to how it is used in the model (a value for each day of each fish's capture history) and then plots a historgram and returns the matrix of covariates (input) z-scored based on the mean and standard deviation of the vector of covariates that correspond with fish's capture histories. 
  scale_mat_func<-function(mat){
    vec<-integer(0)
    for ( i in 1:nrow(dat)){
      vec<-c(vec,mat[(dat[i,"date"]-start_day+1): 
                       ifelse(is.na(dat$TT_MS[i]),
                              end_bon-start_day+1,(dat$TT_MS[i]+dat[i,"date"]-start_day+1)),
                     dat[i,"year"]-2009])
    }
  #  hist(mat)
    (mat-mean(vec))/sd(vec)
    
  }
  pinn_out<-list(log_pinn_Mat=log_pinn_Mat,log_bon_pinn_Mat=log_bon_pinn_Mat)
  
   #use function to scale (i.e. Z-score) (loess predicted) log pinniped count matrix
  log_pinn_Mat_s <- scale_mat_func(log_pinn_Mat)
  
  #scale Bonneville pinniped counts  
  log_bon_pinn_Mat_s<-scale_mat_func(log_bon_pinn_Mat)

  #averages of SCALED Astoria and Bonneville pinnnipeds
  avePin<-scale_mat_func((exp(log_bon_pinn_Mat)+exp(log_pinn_Mat)))#/2
  log_ave_pin<-scale_mat_func(log((exp(log_bon_pinn_Mat)+exp(log_pinn_Mat))))#/2
  
  #load data on Bonneville detection dates for fish tagged as juveniles (datset #2)
  intFile2<-read.csv(here("data","intFile2.csv"))
  #trim to years of Astoria mark-recapture study
  intFile3<-subset(intFile2,detectionYear>=2010 & 
                     detectionYear<=2015)
  rm(intFile2) #delete (for memory conservation)
  
  #drop a few populations that have very few detections
  intFile3<-droplevels(intFile3[!is.na(match(intFile3$Pop,
                                             names(sort(table((intFile3$Pop)))
                                                   [-1:-7]))),])
  
  #drop late pops
 # intFile3<-droplevels(intFile3[is.na(match(intFile3$Pop,
 #                                              c("Pahsimeroi River","East Fork South Fork Salmon","Secesh River","Lostine River","Imnaha River","Upper South Fork Salmon"))),])
                                             # 
  
  #sample sizes of dataset #2
  table(intFile3$Pop,intFile3$detectionYear)
  
  #Because each fish is seen a maximum of one time after it is released, 
  #we can summarize the capture histories with a vectors of release DOY 
  #and the Bonneville DOY, or a value to indicate that it was never seen again. 
  #We already have the vector of Astoria departure date "Astoria_days", 
  #so now I will generate a vector with the Bonneville Arrival days. 
  #For fish that were never seen, I will give them a 
  #Bonneville Arrival day of the day after the last day of the study, 
  #as an indicator.

  bon_Days_Unk_pops<-dat$date+dat$TT_MS-start_day  #"date" is tag date and "TT.Bon" is travel time, so adding them is the Bonneville arrival date
  bon_Days_Unk_pops[dat$surv_ms==0]<-(end_bon-start_day+2) #indicator for fish never seen again
  bon_Days_Unk_pops[dat$surv_ms==1&is.na(dat$TT_MS)]<-(end_bon-start_day+3) # indicator for fish seen only ABOVE Bonneville Dam
  
  #For scaling: vector of all the departure days used in the HMM model (one for each day of each fish's capture history)
  tagDayVec<-rep(dat$date,ifelse(is.na(dat$TT_MS),(end_bon)-(dat$date)+1,dat$TT_MS+1))
 # hist(tagDayVec)
  mean(tagDayVec)
  #*will use log-normalized version of tag day-of-year in model
 # hist(log(tagDayVec))
  length(tagDayVec)
    ############################################################################################
  ##   Reading in environmental covariate values
  ############################################################################################
  ##The data files that I have are weird in that they load in as character data.I have to use sapply to convert them to numeric column by column, Which throws an error when there are any NA's but all the NA's are at the end of each year which I am not using so it is ok.
  
  #read in river temperature at Bonneville data and turn it into numeric column by column (annoying but necessary thing)
  tempDat<-sapply(read.csv(here("data","tempWarr.csv")),
                    function(x)as.numeric(as.character(x)))
  
  #trim to start at starting day of model (not important nwo that I am scaling covariates correctly but I'll leave it int here for now out of lazyness)
  tempDat<-tempDat[start_day:282,7:2]
  
  #make array (necessary for using the "approx" functon for some odd reason i don't undersand)
  tempDat_array<-array(tempDat,dim=c(dim(tempDat),1))
  #linearly interpolate NAs
  tempDat_array[which(is.na(tempDat_array))]<-
    approx(tempDat_array,xout=which(is.na(tempDat_array)))$y
  
  #scale
  tempDat_array<-scale_mat_func(tempDat_array[,,1])
  tempDat_array<-cbind(tempDat_array,rowMeans(tempDat_array),rowMeans(tempDat_array)) #last two are averages, for calculating pinniped effects with other covariates at daily average values
  
  #read in *Discharge* data and turn it into numeric column by column (annoying but necessary thing)
  outflowDat<-sapply( read.csv(here("data","outflowWarr.csv")), 
                      function(x)as.numeric(as.character(x))) #warnings are ok
  #trim to start at starting day of model
  outflowDat<-outflowDat[start_day:282,7:2]
  
  # **I will use a log-transformed version of this covariate because it is otherwise skewed** ##
  logOutflowDat<-log(outflowDat)
  #scale
  outflowDat<- scale_mat_func(outflowDat)
  
  outflowDat<-cbind(outflowDat,rowMeans(outflowDat),rowMeans(outflowDat))#last two are averages, for calculating pinniped effects with other covariates at daily average values
  
  logOutflowDat<- scale_mat_func(logOutflowDat)
  
  
  
  #read in *spill* data
  spillDat<-sapply( read.csv(here("data","spill98-15War.csv"))[,-1],
                    function(x)as.numeric(as.character(x))) #warnings are ok
  spillDat<-spillDat[start_day:282,7:2]
  
  spillDat<-scale_mat_func(spillDat)
  
  spillDat<-cbind(spillDat,rowMeans(spillDat),rowMeans(spillDat))#last two are averages, for calculating pinniped effects with other covariates at daily average values
  

  # End of reading in environmental covariates 
  #########################################################################
  
  ###################################################
  #         End of data loading and wrangling
  ####################################################
  
  ###################################################
  #         Start of running model
  ####################################################
  
  
  #make a list of the data for fitting the model
  data<-list(fit_HMM=1,                #"flag" to fit the HMM
             fit_pop_Ast=1,            #"flag" to fit the population- and year-specific Astoria departure distributions to dataset #2
             penalized_complexity=1,   #"flag" for using penalized complexity priors
             Na=(end_day-start_day+1),                   # number of Astoria departure days to consider
             start_day=start_day,
             Nb=(end_bon-start_day+1),                   # number of Bonneville arrival days to consider
             Nyears=6,                 # number of years that fish were released at Astoria
             
             relDOY=dat$date -start_day,      # tag (Atoria departure) days, but sugtracting 1 so that I can use them as an index in TMB, which starts indexing at 0 not 1.
             
             rel_Year=dat$year-2010,   # year index for each fish released. Subtracting 1 for TMB indexing (starts at 0)
             
             pinCounts=cbind(log_pinn_Mat_s,rowMeans(log_pinn_Mat_s[,1:3]),rowMeans(log_pinn_Mat_s[,4:6])),   #matrix of average log sea lion counts, last two are average of 2010-2012 and 2013-2015
                         phi_x=array(tempDat_array,dim=c(dim(outflowDat),1)),# covariate values for survival 
             
             dayCov=1,                 #is tag-day a covariate in travel time model? yes=1, no=0
             
             X=array(c(  
                       outflowDat,
                       spillDat),dim=c(dim(outflowDat),2)),#other covariate matrices for travel-time 
             
 
             dayScales=c(mean(log(tagDayVec)),
                         sd(log(tagDayVec))),#scaling mean and sd for log-tag day-of-year.
             
             Npops=length(unique(intFile3$Pop)), #number of populations in the data on Bonneville detection dates for fish from known populations
             
             bonDOY_unk_pop=
               bon_Days_Unk_pops-1,   #fate of fish tagged in Astoria: either day of detection at Bonneville or an indicator (value of day after last day of study) that they were never seen again. Subtracting 1 for TMB indexing
             
             relOrigin=dat$clip,      #Origin (hatchery/natural) of fish released in Astoria, based on presence/absence of adipose fin clip.
             
             bonDOY_known_pop=
               intFile3$detectionJulian-start_day,#observation dates of fish (of known population) pasing Bonneville Dam. Subtracting 1 for TMB indexing
             
             pop=as.numeric(as.factor(intFile3$Pop))-1, # index of population for each fish detected at Bonneville Dam, which are used to fit population-specific Astoria departure timing. Subtracting 1 for TMB indexing
             
             bon_Years=intFile3$detectionYear-2009-1, #index of year for each fish detected at Bonneville Dam, which are used to fit population-specific Astoria departure timing. Subtracting 1 for TMB indexing
             exp_rate=c(1,1),              # rate paramaters for exponential priors on standard deviations of normal penalized complexity priors on coefficients. First is for survial model and second is for travel time model. 
             
             ad_mort=ad_mort,
             ad_mig=ad_mig,
             ad_haz=0
  )
  
  
  #Make a list of starting values for parameters
  parameters<-list(log_sigma=runif(6,-1.5,0.5),         #penalized-complexity prior log-standard deviations
                   log_lambda=runif(1,-4,-3),  #Weibull log-sclae
                   log_alpha=runif(1,1,2),          #Weibull log-shape
                   B=rnorm(3,c(1,-.5,1.3),.5),                  #coeffcients for transition (psi) Weibull
                   phi_intercept=rnorm(1,4.8,.5),             #intercept for logit-survival (phi) 
                   phi_betas=rnorm(3,c(-1.2,-.3,-.2),.2),          #coeficcients for logit-survival
                   A_mu=rnorm(length(unique(intFile3$Pop)),1.5,.5),       #population-specific log-mean Astoria departure days
                   A_sigma= rnorm(length(unique(intFile3$Pop)), 1.5,.5), #population-specific log-standard deviation Astoria departure days
                   A_year_beta=rnorm(length(unique(intFile3$detectionYear))-1,0,.1),  # year effects
                   B_year_beta=rnorm(length(unique(intFile3$detectionYear))-1,0,.1))
  
  
  return(list(parameters=parameters,data=data,dat=list(dat,intFile3),pinn=pinn_out))
}

###############             End Data processing         ########################



###########################################
##        Function to plot fits to data                  
###########################################

assessFunc<-function(fit_out, data,start_day,end_day,bon_end, plot_out=F ,print_outs=TRUE,plot_surv=TRUE,plot_TT=TRUE,plot_pop_Ast=FALSE,cex){

  test<-fit_out$mod$report()
  
  dat<-data[[1]]
  intFile3<-data[[2]]
 
  
  sds<-matrix(fit_out$mod_fit$SD$sd[names(fit_out$mod_fit$SD$value)=="all_zeros_prob_or"],ncol=8)
  
####### print likelihood components  
  if(print_outs){
  print(test$l1) #survival and travel time
  print(test$l2)  #population-specific Bonneville detection
  print(test$l3)  # penalized complexity coefficients 
  print(test$l4)  # penalized complexity hyper prior on standard deviations
  print(test$Obj_fun) # total NLL
  }
  ###################Survival by tag date   
  
  if(plot_surv){
  daily_Surv<-tapply(dat$surv_ms, list(dat$date,dat$year,dat$clip), function(x)sum(x)/length(x))
  daily_samp_size<-tapply(dat$surv_ms, list(dat$date,dat$year,dat$clip), function(x)length(x))
  
  all_zeros_prob<-test$all_zeros_prob
  if(plot_out){
    png("surv_plot.png",units="in",height=3.75,width=6,res=300)
  }

   par(mfrow=c(2,4),mar=rep(0,4),oma=c(7,10,2,10),cex=cex)
  #par(cex=1.5,mfrow=c(2,3),mar=c(3,3,3,.5),oma=c(4,4,0,6))
  
  for ( i in c(1:3,7,4:6,8)){
    plot(1,1,type="n",ylim=c(0,1.1),xlim=c(start_day,end_day),
         xlab="",
         ylab="",
         main = "",
         yaxt="n",
         xaxt="n")
    if(i%in%c(4:6,8))axis(1,at=c(91,121,152,182),labels = c("1-Apr","","1-Jun",""),cex.axis=1.5,padj=0)
    if(i%in%c(1,4)) axis(2,cex.axis=1.5)
    if(i<7){
    points(as.numeric(row.names(daily_Surv)),
           1-daily_Surv[,i,1],
           cex=ceiling(daily_samp_size[,i,1]/10),
           col=rgb(0,0,.5,.5),pch=19)
    
    points(as.numeric(row.names(daily_Surv)),
           1-daily_Surv[,i,2],
           cex=ceiling(daily_samp_size[,i,2]/10),
           col=rgb(.5,0,0,.5),pch=19)   
       
    points(start_day:end_day, all_zeros_prob[[2]][,i],type="l",lwd=4,col=rgb(.5,0,0)) 
    }else{
   
      mort_wild<-all_zeros_prob[[1]][,i]
      mort_wild_sd<-sds[,i] 
      polygon(c(start_day:end_day,rev(start_day:end_day)),c(mort_wild+1.96*mort_wild_sd,rev(mort_wild-1.96*mort_wild_sd)),border=FALSE,col=rgb(0,.0,.5,.5))
    }
    
 points(start_day:end_day, all_zeros_prob[[1]][,i],type="l",lwd=4,col=rgb(0,0,.5) )  

 lab<-ifelse(i<7,i+2009,ifelse(i==7,"2010-2012","2013-2015"))
    text(x=(end_day-start_day)*1+start_day ,y=.95,labels = lab, pos=2,cex=1.5,font=2)
    
  }
  points(rep(250,7)-15,c(-.10,.1,.37,.67,1.1,1.53,1.96)/1,
         pch=19,xpd=NA,col=c(rep(rgb(.1,.1,.1,.35),5),rgb(.5,0,0,.5),rgb(0,0,.5,.5)),
         cex=c(1:5,5,5))
  
  text(rep(260,5)-15,c(-.10,.1,.37,.67,1.1,1.53,1.96)/1,
       c("1-10","11-20","21-30","31-40","41-50","Hatchery","Natural"),xpd=NA,pos=4,cex=1.5)
  mtext("Astoria departure day of year",1,4,outer=T,cex=1.2)
  mtext("Mortality",2,4,outer=T,cex=1.2)
  if(plot_out){
   dev.off()
  }
  }
  ###################### Travel Time ################# 

  if(plot_TT){
  tt_calc_func<-function(){
       if(plot_out){
        png("TT_plot.png",units="in",height=4,width=5,res=300)
       }
    par(cex=1.5,mfrow=c(2,3),mar=rep(0,4),oma=c(6,6,3,3))
    for (i in 1:6){
      BA<-test$BcondA[[1]][[i]]/rowSums(test$BcondA[[1]][[i]])
      
      med<-apply(BA,1,function(x) which.min(abs(cumsum(x)-.5)))-(1:(end_day-start_day+1))
      low<-apply(BA,1,function(x)which.min(abs(cumsum(x)-.975)))-(1:(end_day-start_day+1))
      high<-apply(BA,1,function(x)which.min(abs(cumsum(x)-.025)))-(1:(end_day-start_day+1))
   
      
      plot(dat$date[dat$year==2009+i],
           dat$TT_MS[dat$year==2009+i],type="n",xlim=c(start_day,end_day),#range(dat$date),
           ylim=c(0,max(dat$TT_MS,na.rm = T)),main="",
           xlab="",ylab="",pch=19,col=rgb(.3,.3,.3,.6),yaxt="n",xaxt="n")
      if(i%in%c(1,4))(axis(2,at=seq(0,80,by=20)))
      if(i%in%c(4:6))(axis(1,at=c(91,121,152,182),labels = c("1-Apr","","1-Jun","")))
      
      text(x=((end_day-start_day)*0.5)+start_day+20,y=90,pos=1,labels =i+2009,font=2)
     
      polygon(c(start_day:(end_day),(end_day:start_day)),c(low,rev(high)),border=F,col=rgb(.5,0,0,.3))
     points(dat$date[dat$year==2009+i],
                 dat$TT_MS[dat$year==2009+i],type="p",pch=19,col=rgb(.3,.3,.3,.6))
              
              points(start_day:(end_day),med,col=rgb(.5,0,0),type="l",lwd=3)

    }
    mtext("Travel time",2,4,T,xpd=NA,cex=1.25) 
    mtext("Astoria departure day of year",1,4,T,xpd=NA,cex=1.25)
  }
  
  tt_calc_func()
  if(plot_out){
    dev.off()
  }
  }
   ############ population -specific Bonneville Arrival ######
  ##### average across years

  if(plot_pop_Ast){  
    bp<-test$BP

  
  bpAve<-(apply(simplify2array(bp),c(1,2),mean))
  
  if(plot_out){
    png("BonTime_plot.png",units="in",height=8,width=5,res=300)
  }
  
  par(cex=1.5,mfrow=c(6,3),mar=c(2,2,2,0),
      oma=c(2,2,0,4))
  labs<-sub(" .*$", "", levels(out$dat[[2]]$Pop))
  labs[c(3,4,15,16,17)]<-c("EF salmon","EFSF Salmon","Upper GR",  "Upper Salmon", "Upper SF Salmon" )
  
  for ( i in 1:18){
    popSub<-intFile3$detectionJulian[as.numeric(intFile3$Pop)==i]
    histPop<-(hist(popSub,plot=FALSE))
    max(histPop$density)
    hist(popSub,freq=FALSE,xlim=range(out$data$bonDOY_known_pop+out$data$start_day),
         ylim=c(0,max(histPop$density)*1.15),
         main="",
         xlab="",ylab="",yaxt="n"
    )
    mtext(labs[i],3,outer=F,cex=0.8)
    points (start_day:((bon_end)),bpAve[,i],type="l",lwd=2)
    mtext("Bonneville arrival DOY",1,1,outer=T)
    mtext("Average proportion",2,1,outer=T)
  }
  if(plot_out){
    dev.off()
  }
  }
   
 }

########### function to create map that "maps out" entire paramater vectors

create_Map<-function(par_nums,params){
  
  map<-list()
for ( i in par_nums){
  map[[names(params[i])]]<-factor(rep(NA,length(params[[i]])))
}
  return(map)
  }


########### function to fill starting values in from another model fit. 
fixedParams<-function(mod,pars=1:6){
  parameters2<-out$parameters
  parNames=names(parameters2)[pars]
  for(i in parNames){
    
    parameters2[[i]][which(!is.na( mod$env$map[[i]]))]<-mod$env$last.par.best[names(mod$env$last.par.best)==i]
    
    if(is.null( mod$env$map[[i]])){ parameters2[[i]]<-mod$env$last.par.best[names(mod$env$last.par.best)==i]}
    
    
  }
  parameters2
}


############## plot pinnipeds and population-specific survival ###########

#plot pinniped counts
plot_pop_surv<-function(report=result2,start_day=out$data$start_day,end_day=200){
  par(mfrow=c(3,1),mar=c(1,1,0,0),oma=c(2,4,2,2),xpd=NA)  
  layout(matrix(1:3,nrow=3),widths=c(1,1,1),heights=c(1,1,2))
  
  
  
  
  #read in river temperature at Bonneville data and turn it into numeric column by column (annoying but necessary thing)
  tempDat<-sapply(  read.csv(here("data","tempWarr.csv")),
                    function(x)as.numeric(as.character(x)))
  
  #trim to start at starting day of model (not important nwo that I am scaling covariates correctly but I'll leave it int here for now out of lazyness)
  tempDat<-tempDat[out$data$start_day:282,7:2]
  
  #make array (necessary for using the "approx" functon for some odd reason i don't undersand)
  tempDat_array<-array(tempDat,dim=c(dim(tempDat),1))
  #linearly interpolate NAs
  tempDat_array[which(is.na(tempDat_array))]<-
    approx(tempDat_array,xout=which(is.na(tempDat_array)))$y
  
  #plot temp data
  plot(1,1,type="n",ylim=range(tempDat_array,na.rm = TRUE),xlim=c(2009.5,2015.5)+(150/365),xlab="",ylab="Temperature",xaxt="n",yaxt="n")
  #axis(1,at=(2010:2015)+(150/365),labels=F)
  axis(2, at = seq(5,20, by = 5),labels =F )
  axis(2, at = seq(10,20, by = 10),labels =seq(10,20, by = 10) )
  
  for( i in 1:6){
    points(((out$data$start_day:(out$data$start_day+out$data$Na-1))/365)+i+2009,tempDat_array[1:out$data$Na,i,1],type="l",lwd=2)
  }
  
  
  #plot pinniped data
  
  pinnDat1<-exp(out$pinn$log_pinn_Mat)
  pinnDat2<-exp(out$pinn$log_bon_pinn_Mat)
  
  plot(1,1,type="n",ylim=range(pinnDat1[1:out$data$Na,]),xlim=c(2009.5,2015.5)+(150/365),xlab="",ylab="Sea lions",xaxt="n",yaxt="n")
  #axis(1,at=(2010:2015)+(150/365),labels=F)
  axis(2)#, at = seq(-3,2, by = 1.0),labels = seq(-3,2, by = 1.0))
  legend(x="topleft",legend=c("Astoria","Bonneville"),col=c("black","darkgrey"),lty=1,lwd=2)
  
  
  for( i in 1:6){
    points(((out$data$start_day:(out$data$start_day+out$data$Na-1))/365)+i+2009,pinnDat1[1:out$data$Na,i],type="l",lwd=2)
    
    points(((out$data$start_day:(out$data$start_day+out$data$Na-1))/365)+i+2009,pinnDat2[1:out$data$Na,i],type="l",col="darkgrey",lwd=2) 
    
  }
  

  quant_surv<-array(NA,dim=c(18,6,5))
  for ( i in 1:6){
    for ( j in c(14,11,12)){
      quant_surv[j,i,]<- 1- quantile(report[seq(from =j,to=nrow(report),by=18 ),i],probs=c(.025,.25,.5,.75,.975))
    }
  }
  
  
  plot(1,1,type="n",ylim=c(0,1),xlim=c(.5,6.5),xlab="",ylab="Survival",xaxt="n")
  axis(1,at=1:6,labels=2010:2015)
  
  pops<-c("Tucannon River","Minam River","Pahsimeroi River")
  cols<-c('#1b9e77','#d95f02','#7570b3')
  vec<-c(14,11,12)
  for ( i in 1:3){
    points(1:6-.4+(i*.2),quant_surv[vec[i],,3],col=cols[i],pch=19)
    segments(1:6-.4+(i*.2),quant_surv[vec[i],,1],1:6-.4+(i*.2),quant_surv[vec[i],,5],col=cols[i],lwd=2)
  }
  
  legend("bottomleft",legend=paste(pops,c("(Earliest)","(Intermediate)","(Latest)")),col=cols,pch=19,cex=.87,lwd=2)
  
}

################################################
## print coefficients and plot population-specific timing and percent change in survival
###############################################
#function that returns tables of covariate values and produces a plot of average (across years) population-specific Astoria departure timeing and the percent change in average population-specific survial between 2010-2012 and 2013-2015

coef_plot<-function(report,
                    start_day=70,mort,pop_order=c(18,5,10,NA,14,NA,8,11,2,15,6,NA,13,4,17,NA,1,9,NA,7,12,3,16),cols=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'),plot_which=3,cex=1){
  
  labs<-sub(" .*$", "", levels(out$dat[[2]]$Pop))
  labs[c(3,4,15,16,17)]<-c("EF salmon","EFSF Salmon","Upper GR",  "Upper Salmon", "Upper SF Salmon" )
  labs<-labs[pop_order]
  #label MPG's
  MPGs<-c("Upper Columbia","Lower Snake","Grande Ronde/\nImnaha","South Fork\nSalmon","Middle Fork\nSalmon","Upper Salmon")
  
  
  y_loc<--(diff(c(0,which(is.na(pop_order)),(length(pop_order)+1)))/2)+c(which(is.na(pop_order)),(length(pop_order)+1))
  
  col_vec<-rep(cols,each=)
  
  
  A_mu<-report[names(report)=="A_mu"]
  
  A_sigma<-report[names(report)=="A_sigma"]
  
  A_yr<-report[names(report)=="A_year_beta"]

  B_yr<-report[names(report)=="B_year_beta"]
  

  #beta version
  log_mean_Ast<-exp(A_mu)
  log_sd_Ast<-exp(A_sigma)

  
  
  med<-(qbeta(.5,(log_mean_Ast),(log_sd_Ast)))*out$data$Na+start_day
  lower<-(qbeta(.025,(log_mean_Ast),(log_sd_Ast)))*out$data$Na+start_day
  lowerQuant<-(qbeta(.25,(log_mean_Ast),(log_sd_Ast)))*out$data$Na+start_day
  upper<-(qbeta(.975,(log_mean_Ast),(log_sd_Ast)))*out$data$Na+start_day
  upperQuant<-(qbeta(.75,(log_mean_Ast),(log_sd_Ast)))*out$data$Na+start_day
  
  
  
  
  
  num_plots<-ifelse(plot_which>2,2,1)
  
  par(mfrow=c(1,num_plots),mar=c(1,0,0,0),oma=c(2.8,5,1,5),cex=cex)
  
  polygon_func<-function(){
    break_pts<-c(0,which(is.na(pop_order)),20)
    for ( i in 1:3){
      polygon(x=c(-1,365,365,-1),y=rep(break_pts[((i*2)-1):(i*2)],each=2),col="lightgrey",border=FALSE,xpd=FALSE)
      box()
    }
  }
  if(plot_which%in%c(1,3)){
  plot(x=0,y=0,type="n",xlim=range((c(lower,upper))),ylim=c(0,length(pop_order))+.5,yaxt="n",xaxt="n",xlab="",ylab="")
  mtext("Astoria departure",1,2,cex=1, xpd=NA )
  

  polygon_func()
  
  points(med[pop_order],y=1:length(pop_order),pch=19)
  segments((lower)[pop_order],1:length(pop_order),(upper)[pop_order],1:length(pop_order))
  
  segments((lowerQuant)[pop_order],1:length(pop_order)+.2,(lowerQuant)[pop_order],1:length(pop_order)-.2)
 
  segments((upperQuant)[pop_order],1:length(pop_order)+.2,(upperQuant)[pop_order],1:length(pop_order)-.2)
  
  
  
  
  
  
  

  
  
  axis(1,at=c(91,151),labels=c("1-Apr","1-Jun"),cex.axis=.8,padj=-.5)
  axis(2,at=which(!is.na(pop_order)),labels=F)
  text(x=67,y=1:length(pop_order)+.1,
       pos = 2,labels=labs,xpd=NA,srt=0,cex=.6) #121,,182"May","Jul"
  
  if(plot_which!=3){ 
    text(x=max(upper)+1,y=y_loc,
       pos = 4,labels=MPGs,xpd=NA,srt=0,cex=.6)
    axis(4,at=which(is.na(pop_order)),labels=F)
  }
  }
  
  #bootstrap distribution
  if(plot_which%in%c(2,3)){
  mort_qant<-apply(mort, 1, quantile,probs=c(.025,.25,.5,.75,.975))
  
  lower<-mort_qant[1,]
  upper<-mort_qant[5,]
  lowerQuant<-mort_qant[2,]
  upperQuant<-mort_qant[4,]
  med_perc_change<-mort_qant[3,]
  
  
  
  #par(mfrow=c(1,1))
  plot(x=0,y=0,type="n",xlim=c(0,max(upper)),ylim=c(0,length(pop_order))+.5,yaxt="n",xaxt="n",xlab="",ylab="")
  mtext("Mortality",1,2.2,cex=1, xpd=NA )
  polygon_func()
  points(med_perc_change[pop_order],y=1:length(pop_order),pch=19)
  segments((lower)[pop_order],1:length(pop_order),(upper)[pop_order],1:length(pop_order))
  
  segments(lowerQuant[pop_order],1:length(pop_order)+.2,lowerQuant[pop_order],1:length(pop_order)-.2)
  segments(upperQuant[pop_order],1:length(pop_order)+.2,upperQuant[pop_order],1:length(pop_order)-.2)
  
  axis(1,at=seq(0,.3,by=.1),cex.axis=.8,padj=-.5,labels=T)
abline(v=0,lty=2,lwd=1,xpd=T)
  
  axis(4,at=which(is.na(pop_order)),labels=F)
  text(x=max(c(lower,upper))+.02,y=y_loc,
       pos = 4,labels=MPGs,xpd=NA,srt=0,cex=.6)
  
  axis(2,at=which(!is.na(pop_order)),labels=F)
  text(x=-.02,y=1:length(pop_order)+.1,
       pos = 2,labels=labs,xpd=NA,srt=0,cex=.6) #121,,182"May","Jul"
  
}
  
  }


