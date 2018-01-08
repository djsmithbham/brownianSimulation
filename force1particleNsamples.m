% particle with Brownian forcing

clear all;close all;

% SI
a=1e-6;
mu=1e-3;

Nt=6400;
tMax=1.0;

t=linspace(0,tMax,Nt);
dt=t(2)-t(1);

Nsamp=2000;

xp=zeros(Nt,Nsamp);
yp=zeros(Nt,Nsamp);
zp=zeros(Nt,Nsamp);

% R dx/dt = F
% F = sqrt(R/dt)dW
% dx = sqrt(R/dt)dW * dt = sqrt(dt) dW

R  = 6*pi*mu*a;
kB = 1.38e-23;  % Boltzmann's constant
T  = 310;       % absolute temperature

for nt=1:Nt-1
    for ns=1:Nsamp
        F=sqrt(2*kB*T*R/dt)*randn(3,1);
        U=F/R;
        xp(nt+1,ns)=xp(nt,ns)+U(1)*dt;
        yp(nt+1,ns)=yp(nt,ns)+U(2)*dt;
        zp(nt+1,ns)=zp(nt,ns)+U(3)*dt;
    end
end

xsig=std(xp(Nt,:));
ysig=std(yp(Nt,:));
zsig=std(zp(Nt,:));

% can check convergence of the below with Nt
sqrt(xsig^2+ysig^2+zsig^2)



