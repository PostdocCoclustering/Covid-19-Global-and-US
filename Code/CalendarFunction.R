### Plot calendar function:
plot.calendar <- function(df) {
  
  # the data.frame df should have 2 fields: dates and counts
  require(ggplot2); require(lubridate)
  wom <- function(date) { # week-of-month   
    first <- wday(as.Date(paste(year(date),month(date),1,sep="-")))
    return((mday(date)+(first-2)) %/% 7+1)
  }
  df$month <- month(df$dates)
  df$day   <- mday(df$dates)
  
  rng   <- range(as.Date(df$dates))
  rng   <- as.Date(paste(year(rng),month(rng),1,sep="-"))
  start <- rng[1]
  end   <- rng[2]
  month(end) <- month(end)+1
  day(end)   <- day(end)  -1
  
  cal <- data.frame(dates=seq(start,end,by="day"))
  cal$year  <- year(cal$date)
  cal$month <- month(cal$date)
  cal$cmonth<- month(cal$date,label=T)
  cal$day   <- mday(cal$date)
  cal$cdow  <- wday(cal$date,label=T)
  cal$dow   <- wday(cal$date)
  cal$week  <- wom(cal$date)
  
  cal        <- merge(cal,df[,c("dates","counts")],all.x=T) # counts is number of cases/deaths per day
  
  # cal$counts = factor(cal$counts)
  ggplot(cal, aes(cdow, -week, fill = factor(counts))) +
    geom_tile(colour = "white") +
    facet_grid(~cmonth) +
    labs(x="Week of Month",y="",title="",subtitle="",fill=factor("counts")) + # "Groups" are the col.clusters??
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(panel.grid=element_blank())+
    labs(x="",y="")+
    coord_fixed()
}
