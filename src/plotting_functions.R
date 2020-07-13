# functions for plotting model fits and results

#Mark Sorel 7-12-2020

###########################################
##        Function to plot fits to data                  
###########################################
#takes and inputs:model fit object (fit_outs),list of data (data),first DOY for the model and the last DOY when a fish can depart Astoria (end_day) as well as the last day when a fish can arrive at Bonneville Dam (end_bon), and flags for what to print and plot, and point size of plots. Prints plots and parameter values from fitted model object, and saves plots in plot_out=TRUE. 
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

################################################
## print coefficients and plot population-specific timing and sea-lion-associated mortaltiy
###############################################
#takes as input, a report object from a model, first DOY of model, an array of population-specific mortality estimates from a bootstrap procedure (mort), order of populations for plotting, colors for plotting, and an indicator for whether to plot population Astoria departure timing (plot_which=1), morttality (plot_which=2), or both (plot_which=3), and a point size for plots (cex).

#Returns tables of covariate values and produces a plot of average (across years) population-specific Astoria departure timeing and/or the population-specific mortality associated with incresed sea lions between 2010-2012 and 2013-2015

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


