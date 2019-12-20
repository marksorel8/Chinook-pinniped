
#function to make a TS object out of catches
make_ts<-function(myDates,myCounts,main,plot=FALSE,plotDis=FALSE,dis_dates=NULL,dis_vals=NULL){
  
  ts_func<-function(my_dates,my_x){
    
    drop_leap_day<-function(x){
    x[!(format(x,"%m") == "02" & format(x, "%d") == "29")]
  }
  
  out<-ts(
    my_x[match(drop_leap_day(seq.Date(min(my_dates,na.rm=T),
                                          max(my_dates,na.rm=T),by=1)),my_dates)],start=as.numeric(c(format(min(my_dates,na.rm=T), "%Y"),format(min(my_dates,na.rm=T), "%j"))),frequency=365)

  return(out)
  }
  
 out2<-ts_func(myDates,myCounts)
 dis_ts<-NULL
 if(plotDis==TRUE){
 dis_ts<-ts_func(dis_dates,dis_vals)
 }
   
    if(isTRUE(plot)){
    plot_catch(out2,main,plotDis,dis_ts)
  }
  list(catch=out2,dis=dis_ts)

  
}

