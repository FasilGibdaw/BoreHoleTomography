%% MATLAB script to show the position of transmitters and receivers 
% (Figure 1 in the report)
clear all;clc;
dx=1;%step along the horizontal
dy=1;%step along the vertical
x=0:dx:10;
y=0:dy:10;
nx=length(x)-1;ny=length(y)-1;% grid nx by ny
xmin=min(x);xmax=max(x);% end points of the horizontal
ymin=min(y);ymax=max(y);% end points of the vertical
rangex=(xmax-xmin);rangey=(ymax-ymin);
side1=5;     %Number of receivers along side 1
Rxy1=[(xmin+((1:side1)*rangex)/(side1+1))',ymax*ones(side1,1)];
side2=5;     %Number of receivers along side 2
Rxy2=[xmax*ones(side1,1),(ymin+((1:side2)*rangey)/(side2+1))'];
Rxy=[Rxy1;Rxy2];  %Receiver positions alongh side 1 and 2
side3=5;%Number of transmitters along side 3
Txy1=[xmin*ones(side3,1),(ymin+((1:side3)*rangey)/(side3+1))'];
side4=0;%No transmitters along side 4
Txy2=[(xmin+((1:side4)*rangex)/(side4+1))',ymin*ones(side4,1)];
Txy=[Txy1;Txy2];
b=1;          % initalizing number of rays
N_ray=(side1+side2)*(side3+side4);% number of rays from transmitter-receiver pair
ray=zeros(N_ray,4);% initializing ray matrix
for t=1:(side3+side4)
    for r=1:(side1+side2)
           ray(b,1)=Txy(t,1);% x components of the transmitters position
           ray(b,2)=Txy(t,2);% y components of the transmitters position
           ray(b,3)=Rxy(r,1);% x components of the receivers position
           ray(b,4)=Rxy(r,2);% y components of the receivers position
           b=b+1;
       end
end
%% Graphics
% grid plot
for j=1:length(y)
    plot(x,y(j)*ones(length(x)),'k','Linewidth',1)
    hold on
end
hold on
for k=1:length(x)
    plot(x(k)*ones(length(y)),y,'k','Linewidth',1)
    hold on
end
hold on
% Plot position of the transmitters and receivers
plot(Rxy(:,1),Rxy(:,2),'o','markersize',8,'MarkerFaceColor','r')
hold on
plot(Txy(:,1),Txy(:,2),'rs','markersize',8,'MarkerFaceColor','g')
hold on
% Plot ray path between transmitters and receivers
for i=1:N_ray
   h= plot([ray(i,1),ray(i,3)],[ray(i,2),ray(i,4)]);
    hold on
end
xlabel('Position (arbitrary units)')
ylabel('Position (arbitrary units)')
saveas(h,'fig1.png')