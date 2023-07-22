## R script for Tomography
# by Fasil Tesema, November 28 2013, Dessie/Kombolcha, Ethiopia.
# Bahir Dar University, Bahir Dar, Ethiopia
# rm(list=(ls()))# clear all the variables
#setwd("C:/Users/Fasilo/Desktop/MATLAB and R_programs/")
A=as.matrix(read.table('ima.txt'))# the defined unknown 
cc=dim(A)# get the dimension of the input
nx=cc[1];ny=cc[2]# grids in x and y direction (2D)
dx=1;dy=1# step of the grid
x=seq(0,nx,by=dx);y=seq(0,ny,by=dy) #nx=length(x)-1#ny=length(y)-1
xmin=min(x);xmax=max(x);ymin=min(y);ymax=max(y)
rangex=xmax-xmin;rangey=ymax-ymin
side1=5 #number of receivers on side 1 at the surface of the ground
Rxy1=matrix(c(t((xmin+(1:side1)*rangex)/(side1+1)),ymax*(t(rep(1,side1)))),nrow=side1,ncol=2)
side2=5#number of receivers on side 2 (borehole 2)
Rxy2=matrix(c(xmax*(rep(1,side2)),t((ymin+(1:side2)*rangey)/(side2+1))),nrow=side2,ncol=2)
Rxy=rbind(Rxy1,Rxy2)
side3=5#number of transmitters on side 3 (borehole 1)
Txy=matrix(c(xmin*(t(rep(1,side3))),t((ymin+(1:side3)*rangey)/(side3+1))),nrow=side3,ncol=2)

b=1# initializing the variable b
N_ray=(side1+side2)*side3
ray=matrix(0,N_ray,4)# initializing the ray matrix
#########################################################################################
# constructing end points of all the straight rays on the x,y plane
for(t in 1:side3){for (r in 1:(side1+side2)){ray[b,1]=Txy[t,1] 
ray[b,2]=Txy[t,2];ray[b,3]=Rxy[r,1]; ray[b,4]=Rxy[r,2];b=b+1}}
########################################################################################
# constructing the therory matrix for the forward problem.
pixel=seq(1:(nx*ny))# number of pixel
Nr=3000; #Number of points in the discritization along a single ray
G=matrix(0,N_ray,nx*ny)# theory matrix for forward problem
d=0;dr=0;drx=0;dry=0;dd=0# initializing all the required variables in the loop
for (n in 1:N_ray){
    d[n]=sqrt((ray[n,1]-ray[n,3])^2+(ray[n,2]-ray[n,4])^2);
    dr[n]=d[n]/Nr;
    drx[n]=(ray[n,3]-ray[n,1])/Nr;
    dry[n]=(ray[n,4]-ray[n,2])/Nr;
    xp=ray[n,1]+(1:Nr)*drx[n];
    yp=ray[n,2]+(1:Nr)*dry[n];
    ind_x=1+floor(((xp-xmin)/(xmax-xmin))*(rangex-1));
    ind_y=1+floor(((yp-ymin)/(ymax-ymin))*(rangey-1));
    index=(ind_y-1)*nx+ind_x;# index to determine the pixel
    w=table(index); c=as.numeric(w);# counts how many dr are in a paticular pixel
    v=as.numeric(names(w))
    l_ind=c*dr[n];# the total length in a particular pixel
    G[n,v]=G[n,v]+l_ind;# theory matrix with appropriate weighting coefficients
}
dim(A)=c((nx*ny),1)
E=rnorm(N_ray)
m=G%*%A+E# simulated measurement with error
##############################################################################################
#inverse problem
dx1=0.7;dy1=0.7# discretization step of the grid 
x1=seq(0,32,by=dx1);y1=seq(0,52,by=dy1) #
x1min=min(x1);x1max=max(x1);y1min=min(y1);y1max=max(y1)
rangex1=x1max-x1min;rangey1=y1max-y1min
nx1=length(x1)-1;ny1=length(y1)-1
pixel_1=seq(1:(nx1*ny1))# number of pixel != pixel in forward problem
Nr=3000; #Number of points in the discritization along a single ray
G1=matrix(0,N_ray,nx1*ny1)# theory matrix for inversion
d=0;dr=0;drx=0;dry=0;dd=0# initializing all the required variables in the loop
for (n in 1:N_ray){
  d[n]=sqrt((ray[n,1]-ray[n,3])^2+(ray[n,2]-ray[n,4])^2);
  dr[n]=d[n]/Nr;
  drx[n]=(ray[n,3]-ray[n,1])/Nr;
  dry[n]=(ray[n,4]-ray[n,2])/Nr;
  xp=ray[n,1]+(1:Nr)*drx[n];
  yp=ray[n,2]+(1:Nr)*dry[n];
  ind_x=1+floor(((xp-x1min)/(x1max-x1min))*(rangex1-1));
  ind_y=1+floor(((yp-y1min)/(y1max-y1min))*(rangey1-1));
  index=(ind_y-1)*nx1+ind_x;
  w1=table(index); c1=as.numeric(w1);
  v1=as.numeric(names(w1))
  l_ind1=c1*dr[n];
  G1[n,v1]=G1[n,v1]+l_ind1;
}
##################################################################################
# MAP estimation
pr=0.15*rnorm(nx1*ny1);lambda=var(pr)#Gaussian prior
e=0.2*E;sigma=var(e)# expected measurement error
C=rnorm(2*nx1*ny1);cc=var(C);cc=cc*diag(rep(1,2*nx1*ny1));# variance of the innovation
ss=solve(sigma*diag(length(m)))
source('diff_hor.r');source('diff_ver.r')
AA=diff_hor(nx1,ny1);L1=AA[[1]];B=diff_ver(nx1,ny1);L2=B[[1]]# functions
L=rbind(L1,L2)# horizontal and vertical difference
cc=solve(cc)
ss2=t(G1)%*%ss%*%G1+t(L)%*%cc%*%L
ss3=t(G1)%*%ss
xm=solve(ss2,ss3)
xmap=xm%*%m# MAP estimate
dim(xmap)=c(nx1,ny1)# bcak to 2D
image(xmap)
