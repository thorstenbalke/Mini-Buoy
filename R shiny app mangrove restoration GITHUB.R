library(shiny)
library(zoo) ## load package for rollapply function

# Define UI ----
ui <- fluidPage(options(shiny.maxRequestSize = 70*1024^2),
                
                titlePanel(title=div(img(src="my_image.jpg",height="30%", width="30%")))
                ,
                
                sidebarLayout(
                  sidebarPanel(
                    
                    
                    # textInput("buoylocationrest", h3(" Mini buoy location (site of interest)"), 
                    #  value = "Enter text...") ,
                    
                    fileInput("file1", h3("File input: Site of interest")), 
                    
                    
                    
                    # textInput("buoylocationref", h3(" Mini buoy location (reference site)"), 
                    #  value = "Enter text...") ,
                    
                    
                    fileInput("file2", h3("File input: Reference site "))) 
                  
                  ,
                  
                  mainPanel(
                    
                    # Output: Data file ----
                    plotOutput("plot1"),
                    
                    tableOutput('analysis1'),
                    span(textOutput("Warnings1"), style="color:red"),
                    
                    plotOutput('plot2'),
                    tableOutput ('analysis2'),
                    # textOutput ('Warnings2'),
                    span(textOutput("Warnings2"), style="color:red"), 
                    'Difference of conditions: Site of interest - reference site',
                    tableOutput('comparison'),
                    tableOutput("text_Caption")
                    
                  )
                  
                )
                
)








# Define server logic ----
server <- function(input, output, session) {
  
  
  output$plot1<-renderPlot({
    validate(
      need(input$file1, "Please upload 'site of interest' and 'reference' data sets using the browser buttons on the left. 
           \nFor single deployments, please only upload the 'site of interest' data \nPlease refer to 'The Mini Buoy Handbook' for further details.
           \nThe analysis may take a few minutes")
    )
    
    # this is for the MSR145 acceleromter data output with data and time in two separate11 coloumns and with y=-1 acceleration logger in vertical position and y=0 acceleration logger in horizontal position, date and time need to be in a single coloumn
    data1<- read.csv(input$file1$datapath,dec='.',sep=',',skip=38,header=F,fill=T)
    
    
    
    
    data1<-data1[,1:2] 
    names(data1)<-c('date','yacc')
    
    ## set correct time stamp, time zone needs to be defined, GMT is only a placeholder timezone
    data1$date<-strptime(data1$date,format="%Y-%m-%d %H:%M:%S",tz='GMT') 
    
    #derive sampling rate from time stamps
    int<-as.numeric(paste(difftime(data1$date[11],data1$date[10],units = 'secs')))
    rat<-60/int ## samples per minute
    rate1<-rat/60 # rate
    
    ## add 15 min of zeros at start and end to enable deployment and retreival in flooded condition
    start <- strptime(data1$date[1],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    start5min <- start - as.difftime(31*60, units="secs")
    prezero<-data.frame(date= seq(from=start5min,  to=start-int,by = int),yacc=0)
    end<- strptime(data1$date[nrow(data1)],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    end5min<-end + as.difftime(31*60, units="secs")
    postzero<-data.frame(date= seq(from=end+int,  to=end5min,by = int),yacc=0)
    
    #add start and end zero acceleration buffer
    data1<-rbind(prezero,data1,postzero)
    
    
    ## create 15 minute time steps
    
    #current velocity calibration against 15 minute median y acceleraton for aquaculture pond restoration
    curvel15<-function (x) 1.173+((1.059)*x)
    
    
    ##15min running median of acceleration data
    mediany<-rollapply(data1$yacc,width=rat*15,FUN=median,align="center",fill = 0)
    mediany_df<-data.frame(data1,mediany)
    
    
    #select one value every 15 minutes 
    med15_df<-mediany_df[seq(from = 1,to=length(mediany),by=c(rat*15)),]
    
    ##use 15 minute median of yacc <= -0.5 as evidence for submergence 
    velozero<-ifelse(med15_df$mediany<=(-0.5),1,0) #1= flooded , 0 = not flooded
    velocity<-curvel15(med15_df$mediany)*velozero # multiply with velocity i.e. 0 when not flooded
    
    
    ## minimum low tide period of 9*15mins required to separete flood events
    
    veloz60<-rollapply(velozero,9 , sum,align='center',fill=0)#
    
    floodID<-array(data=0,dim=length(velozero))
    
    q=0
    
    for(ii in 1:length(velozero)-1)
      
    {
      ifelse ((veloz60[ii])==0 && (veloz60[ii+1]>=1 ),q<-q+1, q<-q)
      ifelse ((veloz60[ii]>=1),floodID[ii]<-q,floodID[ii]<-0)
      
    }
    
    
    ### for each flood ID determine first and last inundation time step
    predict<-ifelse(velozero==1,'high','low')
    merdat<- data.frame (floodID,med15_df,velocity,velozero,predict)
    
    for (w in 1:max(floodID)) 
    {
      
      s<-min(which(merdat$floodID==w & merdat$predict=='high'))
      e<-max(which(merdat$floodID==w & merdat$predict=='high'))
      
      ifelse(s==e,na<-'na',
             merdat$velozero[min(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('start')
             
      )
      
      ifelse(s==e,na<-'na',
             
             merdat$velozero[max(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('end')
      )
      
    }
    
    
    #row IDs for start/end of inundation
    startend<-data.frame(start=which(merdat$velozero=='start'),end=which(merdat$velozero=='end'))
    
    ### remove 4*15mins start and end of each inundation from hydrodynamic analysis to remove any arifically high currents due to partially inundated loggers
    startend$start<-startend$start+4
    startend$end<-startend$end-4
    #####
    
    #seperate each tide in first and second half (roughly flood and ebb)
    halftide<-startend$start+round((startend$end-startend$start)/2)
    
    sel1<- merdat[startend[1,1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      sel<-merdat[startend[i,1]:startend[i,2],]   
      sel1<-rbind(sel1,sel)
    }
    
    floodtide1<- merdat[startend[1,1]:halftide[1],]
    
    for (i in 2:nrow(startend))
    {
      floodtide<-merdat[startend[i,1]:halftide[i],]   
      floodtide1<-rbind(floodtide1,floodtide)
    }
    
    
    ebbtide1<- merdat[halftide[1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      ebbtide<-merdat[halftide[i]:startend[i,2],]
      ebbtide1<-rbind(ebbtide1,ebbtide)
    }
    
    plot(sel1$date,sel1$velocity,type='n',xaxt='n',ylim=c(0.1,0.5),pch=16,lwd=2,frame=T,xlab=c('Date'),ylab=c('Current velocity (m/s)'),cex.axis=1.3,cex.lab=1.5)
    axis.POSIXct(side = 1,x=sel1$date, format="%b %d",cex=1.5)
    
    points(floodtide1$date,floodtide1$velocity,col='blue',pch=16,cex=1)
    points(ebbtide1$date,ebbtide1$velocity,col='green',pch=16,cex=1)
    abline(h=0.106,col='darkgray',lwd=2)
    
  })
  
  analysis1react<-reactive({
    validate(
      need(input$file1, "")
    )
    data1<- read.csv(input$file1$datapath,dec='.',sep=',',skip=38,header=F,fill=T)
    
    
    
    library(zoo) ## load package for rollapply function
    
    
    data1<-data1[,1:2] 
    names(data1)<-c('date','yacc')
    
    ## set correct time stamp, time zone needs to be defined, GMT is only a placeholder timezone
    data1$date<-strptime(data1$date,format="%Y-%m-%d %H:%M:%S",tz='GMT') 
    
    #derive sampling rate from time stamps
    int<-as.numeric(paste(difftime(data1$date[11],data1$date[10],units = 'secs')))
    rat<-60/int ## samples per minute
    rate1<-rat/60 # rate
    
    ## add 15 min of zeros at start and end to enable deployment and retreival in flooded condition
    start <- strptime(data1$date[1],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    start5min <- start - as.difftime(31*60, units="secs")
    prezero<-data.frame(date= seq(from=start5min,  to=start-int,by = int),yacc=0)
    end<- strptime(data1$date[nrow(data1)],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    end5min<-end + as.difftime(31*60, units="secs")
    postzero<-data.frame(date= seq(from=end+int,  to=end5min,by = int),yacc=0)
    
    #add start and end zero acceleration buffer
    data1<-rbind(prezero,data1,postzero)
    
    
    ## create 15 minute time steps
    
    #current velocity calibration against 15 minute median y acceleraton for aquaculture pond restoration
    curvel15<-function (x) 1.173+((1.059)*x)
    
    
    ##15min running median of acceleration data
    mediany<-rollapply(data1$yacc,width=rat*15,FUN=median,align="center",fill = 0)
    mediany_df<-data.frame(data1,mediany)
    
    
    #select one value every 15 minutes 
    med15_df<-mediany_df[seq(from = 1,to=length(mediany),by=c(rat*15)),]
    
    ##use 15 minute median of yacc <= -0.5 as evidence for submergence 
    velozero<-ifelse(med15_df$mediany<=(-0.5),1,0) #1= flooded , 0 = not flooded
    velocity<-curvel15(med15_df$mediany)*velozero # multiply with velocity i.e. 0 when not flooded
    
    
    ## minimum low tide period of 9*15mins required to separete flood events
    
    veloz60<-rollapply(velozero,9 , sum,align='center',fill=0)#
    
    floodID<-array(data=0,dim=length(velozero))
    
    q=0
    
    for(ii in 1:length(velozero)-1)
      
    {
      ifelse ((veloz60[ii])==0 && (veloz60[ii+1]>=1 ),q<-q+1, q<-q)
      ifelse ((veloz60[ii]>=1),floodID[ii]<-q,floodID[ii]<-0)
      
    }
    
    
    ### for each flood ID determine first and last inundation time step
    predict<-ifelse(velozero==1,'high','low')
    merdat<- data.frame (floodID,med15_df,velocity,velozero,predict)
    
    for (w in 1:max(floodID)) 
    {
      
      s<-min(which(merdat$floodID==w & merdat$predict=='high'))
      e<-max(which(merdat$floodID==w & merdat$predict=='high'))
      
      ifelse(s==e,na<-'na',
             merdat$velozero[min(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('start')
             
      )
      
      ifelse(s==e,na<-'na',
             
             merdat$velozero[max(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('end')
      )
      
    }
    
    
    #row IDs for start/end of inundation
    startend<-data.frame(start=which(merdat$velozero=='start'),end=which(merdat$velozero=='end'))
    
    ### remove 4*15mins start and end of each inundation from hydrodynamic analysis to remove any arifically high currents due to partially inundated loggers
    startend$start<-startend$start+4
    startend$end<-startend$end-4
    #####
    
    #seperate each tide in first and second half (roughly flood and ebb)
    halftide<-startend$start+round((startend$end-startend$start)/2)
    
    sel1<- merdat[startend[1,1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      sel<-merdat[startend[i,1]:startend[i,2],]   
      sel1<-rbind(sel1,sel)
    }
    
    floodtide1<- merdat[startend[1,1]:halftide[1],]
    
    for (i in 2:nrow(startend))
    {
      floodtide<-merdat[startend[i,1]:halftide[i],]   
      floodtide1<-rbind(floodtide1,floodtide)
    }
    
    
    ebbtide1<- merdat[halftide[1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      ebbtide<-merdat[halftide[i]:startend[i,2],]
      ebbtide1<-rbind(ebbtide1,ebbtide)
    }
    
    custa<-summary(sel1$velocity)
    qu<-quantile(sel1$velocity,prob=c(0.25,0.5,0.75))
    
    ##inundation statistics calculations
    
    # time difference in minutes for each inundation 
    inundations<-merdat[merdat$velozero=='start' | merdat$velozero=='end',]
    
    #duration of each inundation
    inundationdur<-data.frame(start=merdat[merdat$velozero=='start',]$date ,end=merdat[merdat$velozero=='end',]$date,duration=0,floodID=merdat[merdat$velozero=='start',]$floodID )
    IDs<-unique(inundations$floodID)
    
    for (ti in 1:nrow(inundationdur))
    { 
      tt<-IDs[ti]
      inundationdur$duration [ti]<- paste(difftime(time1 = inundations[inundations$floodID==tt,] $date [2],time2= inundations[inundations$floodID==tt,] $date [1],units='min'))
    }
    
    
    ## inundaiton free Windows of Opportunity WoO
    
    
    #vecors of start and end of each inundation
    instar<-inundations[inundations$velozero=='start',]
    inend<-inundations[inundations$velozero=='end',]
    
    #all flood IDs
    uniflood<-unique(inend$floodID)
    
    WoO<-array(dim=length(uniflood-1),data=0)
    
    for (wi in 1:c(length(uniflood)-1))
      
    {
      wf<-uniflood[wi]
      wf1<-uniflood[wi+1]
      woendprev<-0
      wostarnew<-0
      woendprev<-inend[inend$floodID==wf,]
      wostarnew<-instar[instar$floodID==wf1,]
      WoO[wi]<-as.numeric( difftime(time1=wostarnew$date,time2=woendprev$date,units='days'))
      
    }
    
    
    #total time difference of measurement period
    timedif<-as.numeric(abs(difftime(data1$date[1],data1$date[nrow(data1)],units = 'days')))
    
    #average inundation durartino per tide
    averagetidedur<-mean(as.numeric(inundationdur$duration))
    
    #average number of tides per day
    tideperday<-nrow(instar)/timedif
    
    #inundation duration per day based on filtered start and end time of tides
    inunminperday<-sum(as.numeric(inundationdur$duration))/timedif
    ##inundation duration per day based on -0.5 threshold only
    inunminperdaytotal=c(60*(sum(velozero)/4)/timedif)
    
    floodstats<-data.frame(measurement_period_days=timedif, average_tide_duration_min=averagetidedur,inundations_per_day=tideperday,max_WoO_days=max(WoO),inundation_per_day_minutes=inunminperdaytotal,median_velocity_ms=qu[2],percentile75_velocity_ms=qu[3],flood_minus_ebb_velocity_ms=c(median(floodtide1$velocity)-median(ebbtide1$velocity)))
    names(floodstats)<-c('Monitoring period (d)','Avg. high tide duration (min)','Flooding frequ. (F/d)','Max. WoO duration (d)','Avg. flooding duration (min/d)','Median current vel. (m/s)','75 percentile velocity (m/s)','Flood-ebb median velocity (m/s)')
    floodstatsa<-t(floodstats)
    
    data.frame(Indicator=names(floodstats),Value=floodstatsa[,1])
    
    
    
    
    
    
    
    #output
  })
  
  output$analysis1<-renderTable({
    analysis1react()
  }, type='html',striped = TRUE, bordered = TRUE,  
  spacing = 'xs', rownames = FALSE,colnames = TRUE,hover = TRUE
  )  
  
  
  
  
  output$Warnings1<-renderText ({
    ana1<-as.data.frame(analysis1react())
    
    warn0<-'Site of interest'
      
    
    #warning if average inundation >1d  
    warn1<-ifelse(ana1[2,2]>=(24*60),'interpretation: Inundation does not follow tidal pattern!',' Interpretation: Inundation is tidal')
    #warning if monitoring period <15 days
    warn2<-ifelse(ana1[1,2]<=(15),'Deployment should be longer',' Deployment duration OK')
    #median velocity not sufficiently different from 75% velocity
    warn3<-ifelse(abs(ana1[6,2]-ana1[7,2])<=(0.01),'current velocities below detection limit','Currents detected OK')
    
    ## inundation not suitable for mangrove >800 <100
    warn4<-ifelse(ana1[5,2]<100||ana1[5,2]>800,'Inundation duration not suitable for mangroves','Inundation duration generally suitable for mangroves')
    # no WoO>1
    warn5<-ifelse(ana1[4,2]<1,'No WoO detected', 'WoO detected')
    
    #high low currents
    warn6<- ifelse(ana1[7,2]<0.15,'Low current velocities','High current velocities')
    
    text1<-paste(warn0,warn1,warn2,warn3,warn4, warn5, warn6,sep=' - ')
    text1
    
  })
  
  
  
  
  
  
  
  ############# reference site
  
  
  
  
  
  
  
  
  
  
  
  
  
  # site1/site of interest/restoration site output table
  
  output$plot2<-renderPlot({
    
    validate(
      need(input$file1, "Data input missing"),
      need(input$file2, "Data input missing")
    )
    
    data1<- read.csv(input$file2$datapath,dec='.',sep=',',skip=38,header=F,fill=T)
    
    library(zoo) ## load package for rollapply function
    
    
    data1<-data1[,1:2] 
    names(data1)<-c('date','yacc')
    
    ## set correct time stamp, time zone needs to be defined, GMT is only a placeholder timezone
    data1$date<-strptime(data1$date,format="%Y-%m-%d %H:%M:%S",tz='GMT') 
    
    #derive sampling rate from time stamps
    int<-as.numeric(paste(difftime(data1$date[11],data1$date[10],units = 'secs')))
    rat<-60/int ## samples per minute
    rate1<-rat/60 # rate
    
    ## add 15 min of zeros at start and end to enable deployment and retreival in flooded condition
    start <- strptime(data1$date[1],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    start5min <- start - as.difftime(31*60, units="secs")
    prezero<-data.frame(date= seq(from=start5min,  to=start-int,by = int),yacc=0)
    end<- strptime(data1$date[nrow(data1)],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    end5min<-end + as.difftime(31*60, units="secs")
    postzero<-data.frame(date= seq(from=end+int,  to=end5min,by = int),yacc=0)
    
    #add start and end zero acceleration buffer
    data1<-rbind(prezero,data1,postzero)
    
    
    ## create 15 minute time steps
    
    #current velocity calibration against 15 minute median y acceleraton for aquaculture pond restoration
    curvel15<-function (x) 1.173+((1.059)*x)
    
    
    ##15min running median of acceleration data
    mediany<-rollapply(data1$yacc,width=rat*15,FUN=median,align="center",fill = 0)
    mediany_df<-data.frame(data1,mediany)
    
    
    #select one value every 15 minutes 
    med15_df<-mediany_df[seq(from = 1,to=length(mediany),by=c(rat*15)),]
    
    ##use 15 minute median of yacc <= -0.5 as evidence for submergence 
    velozero<-ifelse(med15_df$mediany<=(-0.5),1,0) #1= flooded , 0 = not flooded
    velocity<-curvel15(med15_df$mediany)*velozero # multiply with velocity i.e. 0 when not flooded
    
    
    ## minimum low tide period of 9*15mins required to separete flood events
    
    veloz60<-rollapply(velozero,9 , sum,align='center',fill=0)#
    
    floodID<-array(data=0,dim=length(velozero))
    
    q=0
    
    for(ii in 1:length(velozero)-1)
      
    {
      ifelse ((veloz60[ii])==0 && (veloz60[ii+1]>=1 ),q<-q+1, q<-q)
      ifelse ((veloz60[ii]>=1),floodID[ii]<-q,floodID[ii]<-0)
      
    }
    
    
    ### for each flood ID determine first and last inundation time step
    predict<-ifelse(velozero==1,'high','low')
    merdat<- data.frame (floodID,med15_df,velocity,velozero,predict)
    
    for (w in 1:max(floodID)) 
    {
      
      s<-min(which(merdat$floodID==w & merdat$predict=='high'))
      e<-max(which(merdat$floodID==w & merdat$predict=='high'))
      
      ifelse(s==e,na<-'na',
             merdat$velozero[min(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('start')
             
      )
      
      ifelse(s==e,na<-'na',
             
             merdat$velozero[max(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('end')
      )
      
    }
    
    
    #row IDs for start/end of inundation
    startend<-data.frame(start=which(merdat$velozero=='start'),end=which(merdat$velozero=='end'))
    
    ### remove 4*15mins start and end of each inundation from hydrodynamic analysis to remove any arifically high currents due to partially inundated loggers
    startend$start<-startend$start+4
    startend$end<-startend$end-4
    #####
    
    #seperate each tide in first and second half (roughly flood and ebb)
    halftide<-startend$start+round((startend$end-startend$start)/2)
    
    sel1<- merdat[startend[1,1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      sel<-merdat[startend[i,1]:startend[i,2],]   
      sel1<-rbind(sel1,sel)
    }
    
    floodtide1<- merdat[startend[1,1]:halftide[1],]
    
    for (i in 2:nrow(startend))
    {
      floodtide<-merdat[startend[i,1]:halftide[i],]   
      floodtide1<-rbind(floodtide1,floodtide)
    }
    
    
    ebbtide1<- merdat[halftide[1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      ebbtide<-merdat[halftide[i]:startend[i,2],]
      ebbtide1<-rbind(ebbtide1,ebbtide)
    }
    
    plot(sel1$date,sel1$velocity,type='n',xaxt='n',ylim=c(0.1,0.5),pch=16,lwd=2,frame=T,xlab=c('Date'),ylab=c('Current velocity (m/s)'),cex.axis=1.3,cex.lab=1.5)
    axis.POSIXct(side = 1,x=sel1$date, format="%b %d",cex=1.5)
    
    points(floodtide1$date,floodtide1$velocity,col='blue',pch=16,cex=1)
    points(ebbtide1$date,ebbtide1$velocity,col='green',pch=16,cex=1)
    abline(h=0.106,col='darkgray',lwd=2)
    
    
    
  })
  
  analysis2react<-reactive({
    validate(
      need(input$file2, "")
    )
    data1<- read.csv(input$file2$datapath,dec='.',sep=',',skip=38,header=F,fill=T)
    
    library(zoo) ## load package for rollapply function
    
    
    data1<-data1[,1:2] 
    names(data1)<-c('date','yacc')
    
    ## set correct time stamp, time zone needs to be defined, GMT is only a placeholder timezone
    data1$date<-strptime(data1$date,format="%Y-%m-%d %H:%M:%S",tz='GMT') 
    
    #derive sampling rate from time stamps
    int<-as.numeric(paste(difftime(data1$date[11],data1$date[10],units = 'secs')))
    rat<-60/int ## samples per minute
    rate1<-rat/60 # rate
    
    ## add 15 min of zeros at start and end to enable deployment and retreival in flooded condition
    start <- strptime(data1$date[1],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    start5min <- start - as.difftime(31*60, units="secs")
    prezero<-data.frame(date= seq(from=start5min,  to=start-int,by = int),yacc=0)
    end<- strptime(data1$date[nrow(data1)],format="%Y-%m-%d %H:%M:%S",tz='GMT')
    end5min<-end + as.difftime(31*60, units="secs")
    postzero<-data.frame(date= seq(from=end+int,  to=end5min,by = int),yacc=0)
    
    #add start and end zero acceleration buffer
    data1<-rbind(prezero,data1,postzero)
    
    
    ## create 15 minute time steps
    
    #current velocity calibration against 15 minute median y acceleraton for aquaculture pond restoration
    curvel15<-function (x) 1.173+((1.059)*x)
    
    
    ##15min running median of acceleration data
    mediany<-rollapply(data1$yacc,width=rat*15,FUN=median,align="center",fill = 0)
    mediany_df<-data.frame(data1,mediany)
    
    
    #select one value every 15 minutes 
    med15_df<-mediany_df[seq(from = 1,to=length(mediany),by=c(rat*15)),]
    
    ##use 15 minute median of yacc <= -0.5 as evidence for submergence 
    velozero<-ifelse(med15_df$mediany<=(-0.5),1,0) #1= flooded , 0 = not flooded
    velocity<-curvel15(med15_df$mediany)*velozero # multiply with velocity i.e. 0 when not flooded
    
    
    ## minimum low tide period of 9*15mins required to separete flood events
    
    veloz60<-rollapply(velozero,9 , sum,align='center',fill=0)#
    
    floodID<-array(data=0,dim=length(velozero))
    
    q=0
    
    for(ii in 1:length(velozero)-1)
      
    {
      ifelse ((veloz60[ii])==0 && (veloz60[ii+1]>=1 ),q<-q+1, q<-q)
      ifelse ((veloz60[ii]>=1),floodID[ii]<-q,floodID[ii]<-0)
      
    }
    
    
    ### for each flood ID determine first and last inundation time step
    predict<-ifelse(velozero==1,'high','low')
    merdat<- data.frame (floodID,med15_df,velocity,velozero,predict)
    
    for (w in 1:max(floodID)) 
    {
      
      s<-min(which(merdat$floodID==w & merdat$predict=='high'))
      e<-max(which(merdat$floodID==w & merdat$predict=='high'))
      
      ifelse(s==e,na<-'na',
             merdat$velozero[min(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('start')
             
      )
      
      ifelse(s==e,na<-'na',
             
             merdat$velozero[max(which(merdat$floodID==w & merdat$predict=='high'))]<-paste('end')
      )
      
    }
    
    
    #row IDs for start/end of inundation
    startend<-data.frame(start=which(merdat$velozero=='start'),end=which(merdat$velozero=='end'))
    
    ### remove 4*15mins start and end of each inundation from hydrodynamic analysis to remove any arifically high currents due to partially inundated loggers
    startend$start<-startend$start+4
    startend$end<-startend$end-4
    #####
    
    #seperate each tide in first and second half (roughly flood and ebb)
    halftide<-startend$start+round((startend$end-startend$start)/2)
    
    sel1<- merdat[startend[1,1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      sel<-merdat[startend[i,1]:startend[i,2],]   
      sel1<-rbind(sel1,sel)
    }
    
    floodtide1<- merdat[startend[1,1]:halftide[1],]
    
    for (i in 2:nrow(startend))
    {
      floodtide<-merdat[startend[i,1]:halftide[i],]   
      floodtide1<-rbind(floodtide1,floodtide)
    }
    
    
    ebbtide1<- merdat[halftide[1]:startend[1,2],]
    
    for (i in 2:nrow(startend))
    {
      ebbtide<-merdat[halftide[i]:startend[i,2],]
      ebbtide1<-rbind(ebbtide1,ebbtide)
    }
    
    custa<-summary(sel1$velocity)
    qu<-quantile(sel1$velocity,prob=c(0.25,0.5,0.75))
    
    ##inundation statistics calculations
    
    # time difference in minutes for each inundation 
    inundations<-merdat[merdat$velozero=='start' | merdat$velozero=='end',]
    
    #duration of each inundation
    inundationdur<-data.frame(start=merdat[merdat$velozero=='start',]$date ,end=merdat[merdat$velozero=='end',]$date,duration=0,floodID=merdat[merdat$velozero=='start',]$floodID )
    IDs<-unique(inundations$floodID)
    
    for (ti in 1:nrow(inundationdur))
    { 
      tt<-IDs[ti]
      inundationdur$duration [ti]<- paste(difftime(time1 = inundations[inundations$floodID==tt,] $date [2],time2= inundations[inundations$floodID==tt,] $date [1],units='min'))
    }
    
    
    ## inundaiton free Windows of Opportunity WoO
    
    
    #vecors of start and end of each inundation
    instar<-inundations[inundations$velozero=='start',]
    inend<-inundations[inundations$velozero=='end',]
    
    #all flood IDs
    uniflood<-unique(inend$floodID)
    
    WoO<-array(dim=length(uniflood-1),data=0)
    
    for (wi in 1:c(length(uniflood)-1))
      
    {
      wf<-uniflood[wi]
      wf1<-uniflood[wi+1]
      woendprev<-0
      wostarnew<-0
      woendprev<-inend[inend$floodID==wf,]
      wostarnew<-instar[instar$floodID==wf1,]
      WoO[wi]<-as.numeric( difftime(time1=wostarnew$date,time2=woendprev$date,units='days'))
      
    }
    
    
    #total time difference of measurement period
    timedif<-as.numeric(abs(difftime(data1$date[1],data1$date[nrow(data1)],units = 'days')))
    
    #average inundation durartino per tide
    averagetidedur<-mean(as.numeric(inundationdur$duration))
    
    #average number of tides per day
    tideperday<-nrow(instar)/timedif
    
    #inundation duration per day based on filtered start and end time of tides
    inunminperday<-sum(as.numeric(inundationdur$duration))/timedif
    ##inundation duration per day based on -0.5 threshold only
    inunminperdaytotal=c(60*(sum(velozero)/4)/timedif)
    
    floodstats<-data.frame(measurement_period_days=timedif, average_tide_duration_min=averagetidedur,inundations_per_day=tideperday,max_WoO_days=max(WoO),inundation_per_day_minutes=inunminperdaytotal,median_velocity_ms=qu[2],percentile75_velocity_ms=qu[3],flood_minus_ebb_velocity_ms=c(median(floodtide1$velocity)-median(ebbtide1$velocity)))
    names(floodstats)<-c('Monitoring period (d)','Avg. high tide duration (min)','Flooding frequ. (F/d)','Max. WoO duration (d)','Avg. flooding duration (min/d)','Median current vel. (m/s)','75 percentile velocity (m/s)','Flood-ebb median velocity (m/s)')
    
    
    floodstatsa<-t(floodstats)
    
    data.frame(Indicator=names(floodstats),Value=floodstatsa[,1])
    
    
    #output
  })
  
  output$analysis2<-renderTable({
    data.frame(analysis2react())
  }, striped = TRUE, bordered = TRUE,  
  spacing = 'xs', rownames = FALSE
  )  
  
  output$Warnings2<-renderText ({
    
    ana2<-as.data.frame(analysis2react())
    
    warn0a<-'Reference site'
    warn1a<-ifelse(ana2[2,2]>=(24*60),'interpretation: Inundation does not follow tidal pattern!',' Interpretation: Inundation is tidal')
    #warning if monitoring period <15 days
    warn2a<-ifelse(ana2[1,2]<=(15),'Deployment should be longer',' Deployment duration OK')
    #median velocity not sufficiently different from 75% velocity
    warn3a<-ifelse(abs(ana2[6,2]-ana2[7,2])<=(0.01),'current velocities below detection limit','Currents detected OK')
    
    ## inundation not suitable for mangrove >800 <100
    warn4a<-ifelse(ana2[5,2]<100||ana2[5,2]>800,'Inundation duration not suitable for mangroves','Inundation duration generally suitable for mangroves')
    # no WoO>1
    warn5a<-ifelse(ana2[4,2]<1,'No WoO detected', 'WoO detected')
    
    #high low currents
    warn6a<- ifelse(ana2[7,2]<0.15,'Low current velocities','High current velocities')  
    
    text2<-paste(warn0a,warn1a,warn2a,warn3a,warn4a, warn5a, warn6a,sep=' - ')
    text2
    
    
    
  })
  
  
  
  output$comparison<-renderTable({
    
    validate(
      need(input$file1, ""),
      need(input$file2, "No reference site available for comparisson")
    )
    ana1<-as.data.frame(analysis1react())
    ana2<-as.data.frame(analysis2react())
    
    diftab<-data.frame(Indicator=c('Monitoring period (d)','Avg. high tide duration (min)','Flooding frequ. (F/d)','Max. WoO duration (d)','Avg. flooding duration (min/d)','Median current vel. (m/s)','75 percentile velocity (m/s)','Flood-ebb median velocity (m/s)'),Difference=ana1[,2]-ana2[,2])
    diftab
    
    
    
    #difplo<-data.frame(type=as.character(ana1[,1]),difference=diftab)   
    # par(mar=c(10,5,2,0))
    #xx<-barplot(height=difplo[,2],las=2,ylab='Diff. (target - reference)',width=0.5,names.arg = difplo[,1],cex.names = 0.7,xlim=c(0,6))
    # abline(h=0, col='blue',lty=2)
    # text(x=xx, y = 0, label =round( difplo[,2],digits = 2), pos = 1, cex = 0.8, col = "darkblue")      
    
    
  }
  )
  
  
  
  
  
}



# Run the app ----
shinyApp(ui = ui, server = server)