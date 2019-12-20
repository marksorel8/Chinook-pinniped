bon_tot_pin_func<-function(){

make_pinn_ts<-function(dates, counts, plot=FALSE){


library(here)
source(here("src","ts_and_plotting_funcs.R"))


tot_pinn_bonn<-make_ts(dates,counts)$catch

tot_pinn_bonn
}


bonPinn<-read.csv(here("data","Bon PMP window.csv"))
bonPinn$DATE<-as.Date(bonPinn$DATE,format="%m/%d/%Y")

x<-make_pinn_ts(bonPinn$DATE,bonPinn$Total.Pinniped.Abundance,FALSE)

x2<-matrix(c(rep(NA,7),window(x,start=c(2010,8),end=c(2015,365))) ,ncol=6)

x2
}

bon_tot_pin_func()
#abline(v=rep(c(70,210),times=6)+rep((0:5)*365,each=2))
#make_pinn_ts(bonPinn$DATE,bonPinn$ZCA.Abundance,TRUE)
#make_pinn_ts(bonPinn$DATE,bonPinn$EJU.Abundance,TRUE)

