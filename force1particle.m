% particle with Brownian forcing

clear all;close all;

% SI
a=1e-6;
mu=1e-3;

Nt=100;
tMax=1.0;

t=linspace(0,tMax,Nt);
dt=t(2)-t(1);

xp(1)=0;
yp(1)=0;
zp(1)=0;

% R dx/dt = F
% F = sqrt(R/dt)dW
% dx = sqrt(R/dt)dW * dt = sqrt(dt) dW

R  = 6*pi*mu*a;
kB = 1.38e-23;  % Boltzmann's constant
T  = 310;       % absolute temperature

for nt=1:Nt-1
    F=sqrt(2*kB*T*R/dt)*randn(3,1);
    U=F/R;
    xp(nt+1)=xp(nt)+U(1)*dt;
    yp(nt+1)=yp(nt)+U(2)*dt;
    zp(nt+1)=zp(nt)+U(3)*dt;
end

figure(1);clf;plot3(xp,yp,zp);

