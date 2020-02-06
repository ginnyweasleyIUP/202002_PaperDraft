plot_axis <- function(side, at=axTicks(side), ticks=axTicks(side), label="", 
                      lineVal=0.5, lineLab=1.5, adj=0.5, font=1,col="black")
{
  axis(side,at=at,labels=F,col=col)
  mtext(side=side,at=at,ticks,line=lineVal, font=font, col=col)
  mtext(side=side,label,line=lineLab, adj=adj, font=font)
}

DaysSinceToAD <- function (x) 
{
  return(x/360 + 200)
}

# Get solar data
solar <- read.csv("/stacywork/hadcm3/hadcm3_solar_forcing.csv")

# Get volcanic data (4 latitude bands)
path <- "/stacydata/data/squirrel/R_data/CMIP5/Forcing/Crowley2008/"

locs<-c("3090S","030S","030N","3090N")
latsects<-c(-90, -30, 0, 30, 90)
files<-paste("/stacydata/data/squirrel/R_data/CMIP5/Forcing/Crowley2008/ICI5_",locs,"_AOD_c.txt",sep="")

# read the zonal bands of AOD in. The area is equivalent.
D<-matrix(NA,ncol=4,nrow=45001)
for (i in 1:length(files)) {
  tmp<-read.table(files[i])
  D[,i]<-tmp[,2]
}

# time is in the first column
tAD <- list()
tAD$volc<-tmp[,1]
colnames(D)<-locs
rownames(D)<-tAD$volc
# D is now a matrix with the AOD for each band, time in rows (36/year)

volc_mean <- apply(D,1,mean) #D[,1]#[30200:30400]


# PLOTTING ----------------------------------------------------------------
saveToPDF <- F

if(saveToPDF)
{
  filenam <- "044_forcing_data.pdf"
  pdf(file=filenam,width=11,height=6)
}

layout(1,1,1)
par(mar=c(3,3,1,1), oma=c(0,3,0,3), xpd=FALSE)

xlim <- c(800,2010)
yticks <- seq(1365,1366.5,0.5)

#Solar Forcing
solar.x <- DaysSinceToAD(solar$daysSince_1.1.year_hadcm3)
ind.plot <- solar.x>=xlim[1] & solar.x<xlim[2]
plot(solar.x[ind.plot], solar$TSI[ind.plot], xlim=xlim, ylim=c(1364.71,1366.71)+c(0,2.15), type="l", axes=F, xlab="", ylab="",
     panel.first=c(abline(h = yticks, v=axTicks(1), col="grey", lty=3)))
plot_axis(1, label="Time [CE]")
plot_axis(2, at=yticks[c(1,3)], ticks=yticks[c(1,3)], label=expression("TSI [W/m"^"2"*"]"), adj=-0.1)

#Volcanic Forcing
par(new=T)
yticks <- seq(0,0.6,0.2)
ind.plot <- tAD$volc>=xlim[1] & tAD$volc<xlim[2]
plot(tAD$volc[ind.plot],volc_mean[ind.plot], ylim=c(-0.8,0.65), type="l", xlim=xlim, axes=F, xlab="", ylab="",
     panel.first=c(abline(h = yticks, col="grey", lty=3)))
plot_axis(4, at=yticks, ticks=yticks, label="Aerosol Optical Depth", adj=1.1)

box()

if(saveToPDF) dev.off()