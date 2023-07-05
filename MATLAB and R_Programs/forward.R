## R script for Tomography
# by Fasil Tesema, November 28 2013, Dessie/Kombolcha, Ethiopia.
# Phd student at Bahir Dar University, Bahir Dar, Ethiopia
rm(list=(ls()))# clear all the variables
#setwd("C:/Users/Fasilo/Desktop/MATLAB and R_programs/")
nx=10;ny=10# grids in x and y direction (2D)
dx=1;dy=1# step of the grid
x=seq(0,nx,by=dx);y=seq(0,ny,by=dy) #nx=length(x)-1#ny=length(y)-1
xmin=min(x);xmax=max(x);ymin=min(y);ymax=max(y)
rangex=xmax-xmin;rangey=ymax-ymin
side1=5 #number of receivers on side 1 at the surface of the ground
Rxy1=matrix(c(t((xmin+(1:side1)*rangex)/(side1+1)),ymax*(t(rep(1,side1)))),nrow=side1,ncol=2)
side2=5 #number of receivers on side 2 (borehole 2)
Rxy2=matrix(c(xmax*(rep(1,side2)),t((ymin+(1:side2)*rangey)/(side2+1))),nrow=side2,ncol=2)
Rxy=rbind(Rxy1,Rxy2)
side3=5 #number of transmitters on side 3 (borehole 1)
Txy=matrix(c(xmin*(t(rep(1,side3))),t((ymin+(1:side3)*rangey)/(side3+1))),nrow=side3,ncol=2)
##########################################################################
# Graphics
plot(c(Txy[,1],Rxy[,1]),c(Txy[,2],Rxy[,2]),xlab="Position (arbitrary units)",
ylab="Position (arbitrary units)",ylim=c(0,10),pch=NA)
points(Txy[,1],Txy[,2],pch=15,col="green",cex=2)
points(Rxy[,1],Rxy[,2],pch=16,col="red",cex=3)
axis(1,at=0:10)
axis(2,at=0:10)
for (i in 1:length(Txy[,1])){for (j in 1:length(Rxy[,1])){
lines(c(Txy[i,1],Rxy[j,1]),c(Txy[i,2],Rxy[j,2]),type="l",col="blue")
}}
grid(lwd=1,lty="solid",col="black")
box(col="black",lwd=1)
