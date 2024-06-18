####################################################################################
####################################################################################
####################################################################################
####################################################################################
#                        Cases
####################################################################################
library(mgcv)
library(lubridate)
library(cartography)
library(sf)
library(tidyverse)
library(anytime)
library(zoo)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(AICcmodavg)
library(bbmle)
library(MuMIn)
library(patchwork)
library(ggpubr)
library(purrr)
library(plotly)
library(hrbrthemes)
library(ggpubr)
library(gridExtra)
library(grid)
library(lattice)
library(cowplot)
library(plotly)

#============================================================================
#  Figure 1A in main paper  - Shaded time series plot for cases
#============================================================================

ncases<-read.csv("./data/MalawiCovidCases.csv", header = TRUE)

p <- ncases %>% 
  filter(!is.na(ncases$new_cases), !is.na(ncases$dayOfWeek), !is.na(ncases$date)) %>% 
  mutate(date = dmy(date)) 


colfig<-plot_ly(p, x = ~date, y = ~new_cases, type="scatter", mode="markers+text")
fig <- layout(colfig, title = 'COVID-19 cases coloured by variants', xaxis = list(title = 'Time'), 
              yaxis = list(title = 'Number of cases'), 
              shapes = list(
                list(type = "rect",
                     fillcolor = "purple", line = list(color = "purple"), opacity = 0.5,
                     x0 = "2020-04-02", x1 = "2020-12-03", xref = "x",
                     y0 = 0, y1 = 1400, yref = "y"),
                list(type = "rect",
                     fillcolor = "brown", line = list(color = "brown"), opacity = 0.5,
                     x0 = "2020-12-03", x1 = "2021-04-21", xref = "x",
                     y0 = 0, y1 = 1400, yref = "y"),
                list(type = "rect",
                     fillcolor = "green", line = list(color = "green"), opacity = 0.5,
                     x0 = "2021-04-21", x1 = "2021-11-16", xref = "x",
                     y0 = 0, y1 = 1400, yref = "y"),
                list(type = "rect",
                     fillcolor = "orange", line = list(color = "orange"), opacity = 0.5,
                     x0 = "2021-11-16", x1 = "2022-10-19", xref = "x",
                     y0 = 0, y1 = 1400, yref = "y")))
fig



dcases <- data.frame(transform(p, day = as.numeric(date)))
m <- gam(new_cases ~ s(day, bs = "cc", k = 45)+dayOfWeek, data = dcases, method='REML', family = nb(), link=log)

##  Model summary, Table 2 in supplementary material
summary(m)

##  Diagnostic plots, Fig 1 in supplementary material
gam.check(m) 

 


#============================================================================
#  Figure 1C in main paper  -  GAM plot for cases
#============================================================================

p_obj <- plot(m, residuals = T, n = nrow(dcases))
p_obj <- p_obj[[1]] # just one smooth so select the first component
sm_df <- as.data.frame(p_obj[c("x", "se", "fit")])
sm_df <- sm_df %>% add_column(date = dcases$date)

ggplot(sm_df, aes(x = date, y = fit)) +
   geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 1, fill="lightblue") +
  geom_line() + theme_bw() + 
  labs(x = "Date", y = "Partial effect of date")




#============================================================================
#  Figure 2A in main paper  -  Doubling time and growth rate for cases
#============================================================================

DoubleTime1 <- function(dat, timev, aggregate = 14, npts=200, plt=FALSE, subgroup=FALSE, figtitle=""){
  kval <- floor(length(dat)/20)
  kval <- ifelse(kval<10, 10, kval)
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))
  Tv <- timev[1:(length(timev))] 
  DW <- weekdays(as.Date(Tv, origin = "2020-04-02")) #1899-12-30
  datfull <- dat
  dat <- dat[1:(length(dat))]
  
  MGAM <- gam(dat~s(Tv, bs='gp', k=kval)+DW, family=nb())
  
  dow <- rep('Friday', length(timev))  
  newd <- data.frame(Tv=timev, DW=dow)  
  p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
  meanspline <- (exp(mean(c(0,coef(MGAM)[2:7]))+p$fit))
  
  xv<-seq(min(Tv),max(Tv), length=npts)
  dow <- weekdays(as.Date(xv, origin = "2020-04-02"))
  newd <- data.frame(Tv=xv, DW=dow)  
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7  
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  dow <- weekdays(as.Date(xv, origin = "2020-04-02"))
  newd <- data.frame(Tv=xv, DW=dow)  
  X1 <- predict(MGAM, newd,type="lpmatrix")
  Xp <- (X1-X0)/eps  
  Xi <- Xp*0 
  Xi[,1:(kval-1)+7] <- Xp[,1:(kval-1)+7]  
  df <- Xi%*%coef(MGAM)                
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5  
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$time <- as.Date(xv, origin = "2020-04-02")
  if(plt==TRUE){
    xv1<-c(timev, max(timev)+1:14) 
    dow <- weekdays(as.Date(xv1, origin = "2020-04-02"))
    newd <- data.frame(Tv=xv1, DW=dow) # data.frame(Tv=xv, DW=DW)
    p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    dattmp <- data.frame(xval = as.Date(xv1, origin = "2020-04-02"), central = exp(p$fit), 
                         upper= upr, lower=lwr, obs = c(datfull, rep(NA, 14)), meanspl = c(meanspline, rep(NA, 14)))
    
    
    nsim <-1000
    set.seed(42)
    simCIlwr<-apply(array(rnbinom(n=nsim*aggregate, mu=(lwr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE)), c((aggregate),nsim)),2,sum)
    simCIupr<-apply(array(rnbinom(n=nsim*aggregate, mu=(upr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE)), c((aggregate),nsim)),2,sum)
    PredInt <- c(round(sum(exp(p$fit)[length(p$fit)-(aggregate-1):0])), 
                 floor(quantile(simCIlwr, c(0.025), na.rm=T)),ceiling(quantile(simCIupr,0.975, na.rm=T)))             
    resprint <- c(PredInt)
    
  
    
    if(df[npts]-2*df.sd[npts]>0){
      resstring <- paste('Increasing recent trend', sep='')
    }else if(df[npts]+2*df.sd[npts]<0){
      resstring <- paste('Decreasing recent trend', sep='')
    }else{
      if(df[npts]>0){
        resstring <- paste() 
      }else{
        resstring <- paste() 
      }
    }
    
    p1 <- ggplot(data=res, aes(x=time, y=sdt))+
      geom_ribbon(aes(ymin=sdtlow, ymax=sdtup), fill = "orange")+  
      geom_line(color='black') +
      labs(x='Date', title='Derivative of spline arising from GAM', caption=resstring)+
      scale_y_continuous(
        name = expression("Instantaneous growth rate"), 
        sec.axis = sec_axis(~., name = "Doubling Time", 
                            breaks = c(-0.1,-0.05,-0.033,-0.0247, 0,0.0247, 0.033, 0.0495, 0.099, 0.173, 0.3466), 
                            labels=c(-7,-14,-21,-28, 'Infinite',28,21, 14, 7,4,2)), 
        limits = c(min(res$sdtlow), max(res$sdtup)))    
    #########################################################################

    
    print(p1)
    
    res <- list(splinederiv=res, PI=PredInt, splinedaymean=meanspline )
    if(subgroup==TRUE){
      res <- list(splinederiv=res, PI=PredInt, simlow =simCIlwr, simupper =simCIupr, splinedaymean=meanspline )
    }
    
   
  }
  
  res
}

DoubleTime1(dat=dcases$new_cases, timev=seq(1,931,by=1), plt = TRUE)



 

#============================================================================
#  Figure 2c in main paper  -  Coloured growth rate curve for cases
#============================================================================

dt2<-DoubleTime1(dat=dcases$new_cases, timev=seq(1,931,by=1), plt = TRUE)
colfig2<-plot_ly(dt2$splinederiv, x = ~time, y = ~sdt, name = "unemployment", type="scatter", mode="lines+text")

fig2 <- layout(colfig2, title = 'COVID-19 growth rate coloured by variants', xaxis = list(title = 'Time'), 
               yaxis = list(title = 'Gorwth rate'), yaxis2 = list(title = 'Doubling time'),
               shapes = list(
                 list(type = "rect",
                      fillcolor = "purple", line = list(color = "purple"), opacity = 0.5,
                      x0 = "2020-04-02", x1 = "2020-12-03", xref = "x",
                      y0 = -0.12, y1 = 0.2, yref = "y"),
                 list(type = "rect",
                      fillcolor = "brown", line = list(color = "brown"), opacity = 0.5,
                      x0 = "2020-12-03", x1 = "2021-04-21", xref = "x",
                      y0 = -0.12, y1 = 0.2, yref = "y"),
                 list(type = "rect",
                      fillcolor = "green", line = list(color = "green"), opacity = 0.5,
                      x0 = "2021-04-21", x1 = "2021-11-16", xref = "x",
                      y0 = -0.12, y1 = 0.2, yref = "y"),
                 list(type = "rect",
                      fillcolor = "orange", line = list(color = "orange"), opacity = 0.5,
                      x0 = "2021-11-16", x1 = "2022-10-19", xref = "x",
                      y0 = -0.12, y1 = 0.2, yref = "y")))

fig2


















####################################################################################
####################################################################################
####################################################################################
####################################################################################
#                        Deaths
####################################################################################



#============================================================================
#  Figure 1B in main paper  - Shaded time series plot for deaths
#============================================================================

ndeaths<-read.csv("./data/MalawiCovidDeaths.csv", header = TRUE)
 
q<-ndeaths %>% filter(!is.na(ndeaths$new_deaths)) %>% 
  mutate(date =mdy(date)) 
 

colfig3<-plot_ly(q, x = ~date, y = ~new_deaths, type="scatter", mode="markers+text")
fig3 <- layout(colfig3, title = 'COVID-19 deaths coloured by variants ', xaxis = list(title = 'Time'), 
              yaxis = list(title = 'Number of deaths'),
              shapes = list(
                list(type = "rect",
                     fillcolor = "purple", line = list(color = "purple"), opacity = 0.5,
                     x0 = "2020-04-02", x1 = "2020-12-03", xref = "x",
                     y0 = 0, y1 = 74, yref = "y"),
                list(type = "rect",
                     fillcolor = "brown", line = list(color = "brown"), opacity = 0.5,
                     x0 = "2020-12-03", x1 = "2021-04-21", xref = "x",
                     y0 = 0, y1 = 74, yref = "y"),
                list(type = "rect",
                     fillcolor = "green", line = list(color = "green"), opacity = 0.5,
                     x0 = "2021-04-21", x1 = "2021-11-16", xref = "x",
                     y0 = 0, y1 = 74, yref = "y"),
                list(type = "rect",
                     fillcolor = "orange", line = list(color = "orange"), opacity = 0.5,
                     x0 = "2021-11-16", x1 = "2022-10-19", xref = "x",
                     y0 = 0, y1 = 74, yref = "y")))
fig3


ddeaths <- data.frame(transform(q, day = as.numeric(date)))
m1 <- gam(new_deaths ~ s(day, bs = "cc", k = 45)+dayOfWeek, data = ddeaths, method='REML', family = nb(), link=log)
 
##  Model summary, Table 3 in supplementary material
summary(m1)

##  Diagnostic plots, Fig 2 in supplementary material
gam.check(m1) 




#============================================================================
#  Figure 1D in main paper  -  GAM plot for deaths
#============================================================================

q_obj <- plot(m1, residuals = T, n = nrow(ddeaths))
q_obj <- q_obj[[1]]  
sm1_df <- as.data.frame(q_obj[c("x", "se", "fit")])
sm1_df <- sm1_df %>% add_column(date = ddeaths$date)
 
ggplot(sm1_df, aes(x = date, y = fit)) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 1, fill="lightblue") +
  geom_line() + theme_bw() + 
  labs(x = "Date", y = "Partial effect of date")








#============================================================================
#  Figure 2B in main paper  -  Doubling time and growth rate for deaths
#============================================================================

DoubleTime2 <- function(dat, timev, aggregate = 14, npts=200, plt=FALSE, subgroup=FALSE, figtitle=""){
  kval <- floor(length(dat)/20)
  kval <- ifelse(kval<10, 10, kval)
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))
  Tv <- timev[1:(length(timev))] 
  DW <- weekdays(as.Date(Tv, origin = "2020-04-02"))  
  datfull <- dat
  dat <- dat[1:(length(dat))]
  
  MGAM <- gam(dat~s(Tv, bs='gp', k=kval)+DW, family=nb())
  
  dow <- rep('Friday', length(timev))  
  newd <- data.frame(Tv=timev, DW=dow)  
  p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
  meanspline <- (exp(mean(c(0,coef(MGAM)[2:7]))+p$fit))
  
  xv<-seq(min(Tv),max(Tv), length=npts)
  dow <- weekdays(as.Date(xv, origin = "2020-04-02"))
  newd <- data.frame(Tv=xv, DW=dow)  
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7  
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  dow <- weekdays(as.Date(xv, origin = "2020-04-02"))
  newd <- data.frame(Tv=xv, DW=dow)  
  X1 <- predict(MGAM, newd,type="lpmatrix")
  Xp <- (X1-X0)/eps  
  #  print(head(Xp))
  Xi <- Xp*0 
  Xi[,1:(kval-1)+7] <- Xp[,1:(kval-1)+7]  
  df <- Xi%*%coef(MGAM)               
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5  
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$time <- as.Date(xv, origin = "2020-04-02")
  if(plt==TRUE){
    xv1<-c(timev, max(timev)+1:14)#timev #seq(min(Tv),max(Tv), length=npts)+eps
    dow <- weekdays(as.Date(xv1, origin = "2020-04-02"))
    newd <- data.frame(Tv=xv1, DW=dow) # data.frame(Tv=xv, DW=DW)
    p <- predict(MGAM, newd, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    dattmp <- data.frame(xval = as.Date(xv1, origin = "2020-04-02"), central = exp(p$fit), 
                         upper= upr, lower=lwr, obs = c(datfull, rep(NA, 14)), meanspl = c(meanspline, rep(NA, 14)))
    
    
    nsim <-1000
    set.seed(42)
    simCIlwr<-apply(array(rnbinom(n=nsim*aggregate, mu=(lwr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE)), c((aggregate),nsim)),2,sum)
    simCIupr<-apply(array(rnbinom(n=nsim*aggregate, mu=(upr[length(p$fit)-(aggregate-1):0]), size=MGAM$family$getTheta(TRUE)), c((aggregate),nsim)),2,sum)
    PredInt <- c(round(sum(exp(p$fit)[length(p$fit)-(aggregate-1):0])), 
                 floor(quantile(simCIlwr, c(0.025), na.rm=T)),ceiling(quantile(simCIupr,0.975, na.rm=T)))             
    resprint <- c(PredInt)
    
    
    
    if(df[npts]-2*df.sd[npts]>0){
      resstring <- paste('Increasing recent trend', sep='')
    }else if(df[npts]+2*df.sd[npts]<0){
      resstring <- paste('Decreasing recent trend', sep='')
    }else{
      if(df[npts]>0){
        resstring <- paste()#'Plateauing recent trend, \n may be increasing with probability ',round(1-pnorm(0, df[npts],df.sd[npts]), 3), sep='')
      }else{
        resstring <- paste()#'Plateauing recent trend, \n may be decreasing with probability ',round(pnorm(0, df[npts],df.sd[npts]),3), sep='')
      }
    }
    
    p1 <- ggplot(data=res, aes(x=time, y=sdt))+
      geom_ribbon(aes(ymin=sdtlow, ymax=sdtup), fill = "purple")+#+ 
      #      annotate("rect", xmin = as.Date('01/01/2020'), xmax = as.Date('01/01/2021'), ymin = -1, ymax = 1,
      #              alpha = .1,fill = "purple")
      #  geom_ribbon(aes(ymin=-0.1, ymax=0.2, xmin=, xmax=), fill = "purple")+#+ 
      #scale_fill_manual(name = "Dir", values = c("up" = "green", "down" = "red"))+
      geom_line(color='black') +
      labs(x='Date', title='Derivative of spline arising from GAM', caption=resstring)+
      scale_y_continuous(
        name = expression("Instantaneous growth rate"), 
        sec.axis = sec_axis(~., name = "Doubling Time", 
                            breaks = c(-0.1,-0.05,-0.033,-0.0247, 0,0.0247, 0.033, 0.0495, 0.099, 0.173, 0.3466), 
                            labels=c(-7,-14,-21,-28, 'Infinite',28,21, 14, 7,4,2)), 
        limits = c(min(res$sdtlow), max(res$sdtup)))   # theme(
    #########################################################################
    
    
    
    print(p1)
    
    res <- list(splinederiv=res, PI=PredInt, splinedaymean=meanspline )
    if(subgroup==TRUE){
      res <- list(splinederiv=res, PI=PredInt, simlow =simCIlwr, simupper =simCIupr, splinedaymean=meanspline )
    }
    
    
  }
  
  res
}
DoubleTime2(dat=ddeaths$new_deaths, timev=seq(1,926,by=1), plt = TRUE)

 




#============================================================================
#  Figure 2D in main paper  -  Coloured growth rate curve for deaths
#============================================================================

library(plotly)
dt3<-DoubleTime2(dat=ddeaths$new_deaths, timev=seq(1,926,by=1), plt = TRUE)
colfig4<-plot_ly(dt3$splinederiv, x = ~time, y = ~sdt, name = "unemployment", type="scatter", mode="lines+text")
fig4 <- layout(colfig4, title = 'COVID-19 death rate coloured by variants', xaxis = list(title = 'Time'), 
               yaxis = list(title = 'Death rate'), yaxis2 = list(title = 'Doubling time'),
               shapes = list(
                 list(type = "rect",
                      fillcolor = "purple", line = list(color = "purple"), opacity = 0.5,
                      x0 = "2020-04-02", x1 = "2020-12-03", xref = "x",
                      y0 = -0.1, y1 = 0.18, yref = "y"),
                 list(type = "rect",
                      fillcolor = "brown", line = list(color = "brown"), opacity = 0.5,
                      x0 = "2020-12-03", x1 = "2021-04-21", xref = "x",
                      y0 = -0.1, y1 = 0.18, yref = "y"),
                 list(type = "rect",
                      fillcolor = "green", line = list(color = "green"), opacity = 0.5,
                      x0 = "2021-04-21", x1 = "2021-11-16", xref = "x",
                      y0 = -0.1, y1 = 0.18, yref = "y"),
                 list(type = "rect",
                      fillcolor = "orange", line = list(color = "orange"), opacity = 0.5,
                      x0 = "2021-11-16", x1 = "2022-10-19", xref = "x",
                      y0 = -0.1, y1 = 0.18, yref = "y")))

fig4












#============================================================================
#  Figure 3 in main paper  -  Incidence maps for cases
#============================================================================

library(cartography)
library(sf)
library(tidyverse)
library(tmap)

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }
  
  return(labs)
}

 
tcase<-read.csv("./data/GeographicData.csv", header = TRUE)
gdata=st_read("./shps/sdr_subnational_boundaries2.shp")
mapdata<-left_join(gdata,tcase,by="DHSREGSP")  
mapdata1<-mapdata %>% filter(!is.na(mapdata$TotalCases))



#============================================================================
#  Figure 3A in main paper  -  Total reported cases
#============================================================================

brks <- c(0, 500, 1000, 2000, 3000, 10000, 20000, 25000)  

labs <- create_labels(brks, greater = F)
pal <- colorRampPalette(c("green", "red"))(7)
map1 <- tm_shape(mapdata1) +
  tm_polygons(col = "TotalCases", palette = pal, labels = labs, breaks = brks, style = "fixed", legend.show = F) +
  tm_add_legend(type = "fill",
                labels = labs, col = pal,
                title = "Total Reported Cases",
                is.portrait = T) +
  tm_layout(legend.outside = F, inner.margins = c(0.05, 0.25, 0.05, 0.35)) + # if you want the legend to be outside, change the legend.outside=T
  tm_compass(position = c("RIGHT", "BOTTOM")) +
  tm_scale_bar(position = c("LEFT", "BOTTOM"))
map1



 

#============================================================================
#  Figure 3B in main paper  -  Total population
#============================================================================

mapdata2<-mapdata %>% filter(!is.na(mapdata$Population))
brks <- c(0, 100000, 200000, 400000, 1000000, 1200000, 1500000, 3000000)  

labs <- create_labels(brks, greater = F)
pal <- colorRampPalette(c("green", "red"))(7)   
map2 <- tm_shape(mapdata2) +
  tm_polygons(col = "Population", palette = pal, labels = labs, breaks = brks, style = "fixed", legend.show = F) +
  tm_add_legend(type = "fill",
                labels = labs, col = pal,
                title = "Total Population",
                is.portrait = T) +
  tm_layout(legend.outside = F, inner.margins = c(0.05, 0.25, 0.05, 0.35)) +  
  tm_compass(position = c("RIGHT", "BOTTOM")) +
  tm_scale_bar(position = c("LEFT", "BOTTOM"))
map2






#============================================================================
#  Figure 3C in main paper  -  Case incidence
#============================================================================

mapdata3<-mapdata %>% filter(!is.na(mapdata$Incidence))
brks <- c(0,1,2,3,4,5,6,7,10,20)  

labs <- create_labels(brks, greater = F)
pal <- colorRampPalette(c("green", "red"))(9)   
map3 <- tm_shape(mapdata3) +
  tm_polygons(col = "Incidence", palette = pal, labels = labs, breaks = brks, style = "fixed", legend.show = F) +
  tm_add_legend(type = "fill",
                labels = labs, col = pal,
                title = "Case Incidnce per 1,000 people",
                is.portrait = T) +
  tm_layout(legend.outside = F, inner.margins = c(0.05, 0.25, 0.05, 0.35)) + # if you want the legend to be outside, change the legend.outside=T
  tm_compass(position = c("RIGHT", "BOTTOM")) +
  tm_scale_bar(position = c("LEFT", "BOTTOM"))
map3




#============================================================================
#  Figure 3 in supplementary material - CFR
#============================================================================
##### Posterior distributions on the same plot ##################################################


x <- c(0:10000)/200000
y1 <- dbeta(x,185,6044-185)   # Other 
y2 <- dbeta(x,957,27933-957)  # Beta
y3 <- dbeta(x,1160,27872-1160) # Delta
y4 <- dbeta(x,381,26219-381)  # Omicron

dat<-data.frame("CFR" = 100*x, "Other" = y1, "Beta" = y2, "Delta" = y3, "Omicron" = y4)
head(dat)

ggplot(data=dat, aes(x=CFR)) + 
  geom_line(aes(y=Other, color='Other'), linewidth=1.5, alpha=0.8) +
  geom_line(aes(y=Beta, color='Beta'),linewidth=1.5, alpha=0.8) +
  geom_line(aes(y=Delta, color='Delta'),linewidth=1.5, alpha=0.8) + 
  geom_line(aes(y=Omicron, color='Omicron'),linewidth=1.5, alpha=0.8)+
  #ggtitle("Posterior distribution of CFR for the variants")+
  xlab("CFR (%)")+
  ylab("Posterior density")+
  scale_colour_manual("COVID-19 variants", 
                      breaks = c("Other", "Beta", "Delta", "Omicron"),
                      values = c("Other"="Purple", "Beta"="Brown", 
                                 "Delta"="Green", "Omicron"="Orange"))

ggsave('./cfrs.pdf', width = 5, height = 3)


#########################################################################################
#########################################################################################
########################################################################################
##  chi-square test and plot Residuals
#########################################################################################
library(corrplot)
library("gplots")

#VARIANT DATA
variants<-matrix(c(4,0,0,2,180,78,36,198,383,42,5,139,213,21,37,2,30,3,4,45), ncol=4, byrow=TRUE) #with missing data-new data
colnames(variants)<-c("South", "Central", "North", "No region")
rownames(variants)<-c("Alpha", "Beta", "Delta", "Omicron", "Other")
 
chiv<-chisq.test(variants, simulate.p.value=TRUE,correct=FALSE)
chiv


#============================================================================
#  Figure 4A in supplementary material - contingency plot
#============================================================================

dt <- as.table(as.matrix(variants))
balloonplot(t(dt), main ="", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)


#============================================================================
#  Figure 5F in main paper - Residual plot
#============================================================================

round(chiv$residuals, 3)
corrplot(chiv$residuals, is.cor = FALSE)


#============================================================================
#  Figure 4B in supplementary material - Cell contribution to Chi-square test
#============================================================================

contrib <- 100*chiv$residuals^2/chiv$statistic
round(contrib, 3)
corrplot(contrib, is.cor = FALSE)




#============================================================================
# Table 4 in supplementary material - proportions for variants 
#============================================================================

if(!require(DescTools)){install.packages("DescTools")}
if(!require(PropCIs)){install.packages("PropCIs")}
library(DescTools)
library(PropCIs)
library(car)


#Alpha
prop.test(x=4, n=6, conf.level=0.95)#inappropriate, fails np and nq >5

#Beta
prop.test(x=180, n=492, conf.level=0.95)  
prop.test(x=78, n=492, conf.level=0.95) 
prop.test(x=36, n=492, conf.level=0.95) 
prop.test(x=198, n=492, conf.level=0.95)  

#Delta
prop.test(x=383, n=569, conf.level=0.95) 
prop.test(x=42, n=569, conf.level=0.95) 
prop.test(x=5, n=569, conf.level=0.95)  
prop.test(x=139, n=569, conf.level=0.95) 

#Omicron 
prop.test(x=213, n=273, conf.level=0.95) 
prop.test(x=21, n=273, conf.level=0.95) 
prop.test(x=37, n=273, conf.level=0.95) 
prop.test(x=2, n=273, conf.level=0.95)  

#Other 
prop.test(x=30, n=82, conf.level=0.95) 
prop.test(x=3, n=82, conf.level=0.95)  
prop.test(x=4, n=82, conf.level=0.95)  
prop.test(x=45, n=82, conf.level=0.95)  


























