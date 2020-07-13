# a functions for processing data and creating initial paramater values for Lower Columbia CHinook survival model

#Mark Sorel 7-12-2020

###################################################
#              Data loading and wrangling
####################################################
#takes as input the first DOY for the model and the last DOY when a fish can depart Astoria (end_day) as well as the last day when a fish can arrive at Bonneville Dam (end_bon). Returns a list of initial paramater values for model, data for model, and some data that can be used for plotting. 
dataProcess<-function(start_day=75,end_day=200,end_bon=210){ #function to read and process data, returns list of data for model, and a list of initial parameters for model, as well as data for plotting
  
  library(here)
 
#-----------------------------------------------------------------------------------------------  
#     Load and process survival and travel time data data
  #-----------------------------------------------------------------------------------------------  
  
  #read data with Astoria tagging, survival, and Bonneville detection date (for survivors) data (Dataset #1)
  dat<-read.csv(here("data","surv_dat_MS.csv"))
  
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
  
  
  #-----------------------------------------------------------------------------------------------  
  #     Load and process population-specific Bonneville arrival data 
  #-----------------------------------------------------------------------------------------------  
  
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
  
  
  
  #-----------------------------------------------------------------------------------------------  
  #     Load and process covariate data
  #-----------------------------------------------------------------------------------------------  
  
#Because a covariate value will be used for each day of each fish's capture history, we need to scale the covariates based on how they are used in the model. This function creates a vector of a covariate equivalent to how it is used in the model (a value for each day of each fish's capture history) and then returns the matrix of covariates that have been z-scored based on the mean and standard deviation of the vector of covariates that correspond with all fish's capture histories across all days and years. 
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
  

  #******load and process sea lion data******
  if(file.exists(here("data","sea_lion_loess.csv"))){        # Check if loess smoothed values have been saved
    predPinn<-read.csv(here("data","sea_lion_loess.csv"))[,1]
    }else
    {
      
  #load Astoria sea lion count data
  pinnDat<-read.csv(here("data","pinnCounts.csv"))
  
  #create a time series with NAs for days with no pinniped counts 
  ##create sequence of zeros
  dateSeq<-paste(rep(2010:2015,each=365),rep(1:365,times=length(2010:2015)))
  ##make data frame
  pinnTS<-data.frame(year=rep(2010:2015,each=365),doy=rep(1:365,times=length(2010:2015)),
                     count=pinnDat$Sum.of.Max.of.Count[match(dateSeq,paste(pinnDat$year,pinnDat$julian))])
   
  #fit a loess model to the pinniped counts to interpolate missing values.
  loessTest<-loess(log(pinnTS$count+1)~c(1:length(pinnTS$count)),span=.02)
  #make predictions for daily pinniped counts from loess model
  predPinn<-ts(predict((loessTest),
                       newdata =data.frame( x=c(1:length(pinnTS$count)))),
               start=c(2010,1),frequency = 365)
  write.csv(predPinn,file=here("data","sea_lion_loess.csv"),row.names = FALSE)
    }
  
  #make a matrix of Astoria sea lion counts (with a column for each year) 
  log_pinn_Mat<-matrix((predPinn),ncol=6)[start_day:282,]
 
   pinn_out<-list(log_pinn_Mat=log_pinn_Mat) #list for returning
  
   #use function to scale (i.e. Z-score) (loess predicted) log pinniped count matrix
  log_pinn_Mat_s <- scale_mat_func(log_pinn_Mat)
  
  
  #******Tag DOY for travel time model******
  #For scaling, creates a vector of all the departure days used in the HMM model (one for each day of each fish's capture history)
  tagDayVec<-rep(dat$date,ifelse(is.na(dat$TT_MS),(end_bon)-(dat$date)+1,dat$TT_MS+1))

    ############################################################################################
  ##   Reading in environmental covariate values
  ############################################################################################
  ##The data files that I have are weird in that they load in as character data.I have to use sapply to convert them to numeric column by column.
  
  #read in ****river temperature**** at Bonneville data and turn it into numeric column by column (annoying but necessary thing)
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
  
  #read in *****Discharge***** data and turn it into numeric column by column (annoying but necessary thing)
  outflowDat<-sapply( read.csv(here("data","outflowWarr.csv")), 
                      function(x)as.numeric(as.character(x))) #warnings are ok
  #trim to start at starting day of model
  outflowDat<-outflowDat[start_day:282,7:2]
  #scale
  outflowDat<- scale_mat_func(outflowDat)
  
  outflowDat<-cbind(outflowDat,rowMeans(outflowDat),rowMeans(outflowDat))#last two are averages, for calculating pinniped effects with other covariates at daily average values
  
  
  #read in ******spill******** data
  spillDat<-sapply( read.csv(here("data","spill98-15War.csv"))[,-1],
                    function(x)as.numeric(as.character(x))) #warnings are ok
  spillDat<-spillDat[start_day:282,7:2]
  
  spillDat<-scale_mat_func(spillDat)
  
  spillDat<-cbind(spillDat,rowMeans(spillDat),rowMeans(spillDat))#last two are averages, for calculating pinniped effects with other covariates at daily average values
  

  # End of reading in environmental covariates 
  
  ###################################################
  #         End of data loading and wrangling
  ####################################################
 
  
   #-----------------------------------------------------------------------------------------------  
  #     Bundle data for model and make initial parameter values
  #-----------------------------------------------------------------------------------------------  

  
  
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
             
             X=array(c(outflowDat,
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
             
             ad_mort=0,       # flag for whether to calculate delta method standard errors for survival model
             ad_mig=0,         # flag for whether to calculate delta method standard errors for population-specifica Astoria depature timeing model
             ad_haz=0               # flag for whether to calculate delta methods tandard erros for paramaters in daily Bonnveille passage probabilities (travel time) model
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
}#end function

###############             End Data processing  function        ########################



