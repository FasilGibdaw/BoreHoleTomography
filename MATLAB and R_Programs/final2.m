%% MATLAB script for Borehole Tomography simulation
% Bahir Dar University, Ethiopia, November 22 2013,
% @Fasil Tesema, email: fasgibe@yahoo.com
%% All the inputs for the problem
clear all;clc;
A=double(imread('endi5.png'));% input image
ima=A(:,:,1); ima=-(ima-128);% true model
dx=1;dy=1; %sampling step along x and y axis
ima=ima(1:dx:end,1:dy:end); % sampled model
[ny,nx]=size(ima);% determines the grid along the x and y coordinates
ima=reshape(ima,nx*ny,1);% making the measurement a columon vector
x=0:nx;% grid points along the x axis
y=0:ny;% grid points along the y axis
xmin=min(x);xmax=max(x);% end points of the x axis
ymin=min(y);ymax=max(y);% end points of the y axis
%% Defining receivers and transmitters positions
side1=5;     %Number of receivers along side 1 at the surface of the ground
rangex=xmax-xmin;rangey=ymax-ymin;
Rxy1=[(xmin+((1:side1)*rangex)/(side1+1))',ymax*ones(side1,1)];
side2=5;     %Number of receivers along side 2 (borehole 2)
Rxy2=[xmax*ones(side2,1),(ymin+((1:side2)*rangey)/(side2+1))'];
Rxy=[Rxy1;Rxy2];  %Receiver positions alongh side 1 and 2
side3=5;%Number of transmitters along side 3 (borehole 1)
Txy1=[xmin*ones(side3,1),(ymin+((1:side3)*rangey)/(side3+1))'];
side4=0;
Txy2=[(xmin+((1:side4)*rangex)/(side4+1))',ymin*ones(side4,1)];
Txy=[Txy1;Txy2];
%% constructing rays with corresponding x and y coordinates
b=1;          % initalizing number of rays
N_ray=(side1+side2)*(side3+side4);% total number of rays 
% initializing the ray matrix with coordinates of end points of all the number of rays
ray=zeros(N_ray,4);% here 4 indicates p1(x1,y1) and p2(x2,y2) of the ray
for t=1:(side3+side4) %index for number of transmitters 
    for r=1:(side1+side2) %index for number of receivers
           ray(b,1)=Txy(t,1);% x components of p1
           ray(b,2)=Txy(t,2); % y components of p1
           ray(b,3)=Rxy(r,1);% x components of p2
           ray(b,4)=Rxy(r,2);% y components of p2
           b=b+1;
       end
end
pixel=1:nx*ny; % number of pixel in one vector
%% Calculating the theory matrix for the forward problem to simulate the measurements
Nr=3000; %number of disretization along each ray
G=zeros(N_ray,length(pixel));% initializing the theory matrix
for n=1:N_ray
    d(n)=sqrt((ray(n,1)-ray(n,3))^2+(ray(n,2)-ray(n,4))^2);%length of a ray
    dr(n)=d(n)/Nr;% step length of a ray
    drx(n)=(ray(n,3)-ray(n,1))/Nr;% x coordinate of step length
    dry(n)=(ray(n,4)-ray(n,2))/Nr;% y coordinate of step length
    xp=ray(n,1)+(1:Nr)*drx(n);% number of x coordinates alon a ray
    yp=ray(n,2)+(1:Nr)*dry(n);% number of y coordinates alon a ray
    ind_x=1+floor(((xp-xmin)/(rangex))*(nx-1));
    ind_y=1+floor(((yp-ymin)/(rangey))*(ny-1));
    index=(ind_y-1)*nx+ind_x;
    dd=unique(index);
    cc=hist(index,dd);
    l_ind=cc*dr(n);
    G(n,dd)=G(n,dd)+l_ind;
end
%% Generating the simulated measurement m with error E 
E=13.3*randn(N_ray,1);%assumed error in the measurement
m=G*ima+E; %the measurement
%% Disretization for the inverse problem
dx1=0.7;dy1=0.7;% discretization step different from dx and dy from above
x1=0:dx1:nx;% grid points along the x axis
y1=0:dy1:ny;% grid points along the y axis
nx1=length(x1)-1;ny1=length(y1)-1;
x1min=min(x1);x1max=max(x1);% end points of the x axis
y1min=min(y1);y1max=max(y1);% end points of the y axis
rangex1=x1max-x1min;rangey1=y1max-y1min;
pixel_1=1:nx1*ny1; % number of pixel in one vector
%% Calculating the therory matrix for the inversion problem
Nr=3000; %number of disretization along each ray
G1=zeros(N_ray,length(pixel_1));% initializing the theory matrix
for n=1:N_ray
    d(n)=sqrt((ray(n,1)-ray(n,3))^2+(ray(n,2)-ray(n,4))^2);%length of a ray
    dr(n)=d(n)/Nr;% step length of a ray
    drx(n)=(ray(n,3)-ray(n,1))/Nr;% x coordinate of step length
    dry(n)=(ray(n,4)-ray(n,2))/Nr;% y coordinate of step length
    xp=ray(n,1)+(1:Nr)*drx(n);% number of x coordinates alon a ray
    yp=ray(n,2)+(1:Nr)*dry(n);% number of y coordinates alon a ray
    ind_x=1+floor(((xp-x1min)/(rangex1))*(nx1-1));
    ind_y=1+floor(((yp-y1min)/(rangey1))*(ny1-1));
    index=(ind_y-1)*nx1+ind_x;
    dd=unique(index);
    cc=hist(index,dd);
    l_ind=cc*dr(n);
    G1(n,dd)=G1(n,dd)+l_ind;
end
%% Constructing the prior (difference prior) and MAP estimate
e=0.2*E;sigma=var(e); % assumed measurement error covariance
L1=diff_ver(nx1,ny1);%function to compute difference matrix along horizontal(L1)  
L2=diff_hor(nx1,ny1);%function to compute difference matrix along vertical(L2)
LL=[L2;L1];% difference matrix (horizontal + vertical)
C=1.2*randn(1,2*nx1*ny1);%virtual measurement error (uncertainity in the difference)
cc=var(C);cc=cc*eye(2*nx1*ny1);%covariance of the uncertainity in the difference
Xm=inv((G1'*inv(sigma*eye(length(m)))*G1)+LL'*inv(cc)*LL)*G1'*inv(sigma*eye(length(m)));
XmapN=Xm*m;% MAP estimate
XmapN=reshape(XmapN,ny1,nx1);% back to the orginal dimension (2D)
ima=reshape(ima,ny,nx);% back to the orginal dimension (2D) of assumed true model
%% Graphics
%figure('PaperPosition',rect = [left, bottom, width, height])
x0=0;
y0=0;
width=900;
height=400
set(gcf,'position',[x0,y0,width,height])
subplot(131)
pcolor(ima);title('\bf (a) True model (625 pixels)');
 subplot(132)
pcolor(ima);
 hold on
DN=1;% DN is the step size to select rays to be plotted
NN=1:DN:N_ray;
for i=1:length(NN)
    plot([ray(NN(i),1),ray(NN(i),3)],[ray(NN(i),2),ray(NN(i),4)],'c')
    hold on
end

axis([xmin xmax ymin ymax])
title('\bf (b) True model with rays')
subplot(133)
pcolor(XmapN);
title('\bf (c) MAP estimate (1225 pixels)')



