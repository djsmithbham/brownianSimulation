% elastic force on 2 particles connected with a spring. Hydrodynamic drag
% but no hydrodynamic interaction

% SI
a=1e-6;
mu=1e-3;
k=2e-7; % not sure about this! spring constant
b0=3e-6; % equilibrium bond length
R=6*pi*mu*a;

kB = 1.38e-23;  % Boltzmann's constant
T  = 310;       % absolute temperature

Np=6;
Nt=200;
tMax=2;

t=linspace(0,tMax,Nt);
dt=t(2)-t(1);

xp=zeros(3,Np,Nt);

for np=1:Np
    xp(1,np,1)=4e-6*(np-1);
    xp(2,np,1)=1e-6*(np-1).^2;
end
f=zeros(3,Np);

%-------------------------------------------------------------------------
% elastic potential is   1/2 k (b-b0)^2
%              = 1/2 k * 2 * (b-b0) * db/dxj
% but b = sqrt((xa-xb)^2+(ya-yb)^2+(za-zb)^2)
% so  db/dxj = 1/2 * 2 (xaj-xbj) / b
%            = (xaj-xbj) / b
% so force is - k (b-b0) (xaj-xbj) / b

% for a string of particles, have force from 'left' and force from 'right'
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Bending potential...
% 
% try potential of sum kb (1-cos(th_j)) where j=1...Np-2
%  = kb sum (1 - (xj+2-xj+1).(xj+1-xj) / |xj+2-xj+1| / |xj+1-xj| )
%
% then force on particle 1 is
%
%  = kb (d/dx1) [ (x3-x2).(x2-x1) / |x3-x2| / |x2-x1| ]
%  = kb [ -(x3-x2) / |x3-x2| / |x2-x1| 
%         -(x3-x2).(x2-x1) / |x3-x2| . (x1-x2) / |x2-x1|^3 ]
%
% have checked above differentiation - see checkdiff.m
% 
% force on particle 2 has two contributions, from th_1 and from th_2
%
%  = kb (d/dx2) [ (x3-x2).(x2-x1) / |x3-x2| / |x2-x1| 
%                +(x4-x3).(x3-x2) / |x4-x3| / |x3-x2| ]
%  = kb [ -(x2-x1) / |x3-x2| / |x2-x1|
%         +(x3-x2) / |x3-x2| / |x2-x1|
%         -(x3-x2).(x2-x1) / |x3-x2| . (x2-x1) / |x2-x1|^3
%         -(x3-x2).(x2-x1) / |x2-x1| . (x2-x3) / |x3-x2|^3
%         -(x4-x3) / |x4-x3| / |x3-x2| 
%         -(x4-x3).(x3-x2) / |x4-x3| . (x2-x3) / |x3-x2|^3 ]
%
%
% force on particle 3 has three contributions, from th_1, th_2 and th_3
%
% = kb (d/dx3) [ (x3-x2).(x2-x1) / |x3-x2| / |x2-x1| 
%               +(x4-x3).(x3-x2) / |x4-x3| / |x3-x2|
%               +(x5-x4).(x4-x3) / |x5-x4| / |x4-x3| ]
%
% = kb [  (x2-x1) / |x3-x2| / |x2-x1|
%        -(x3-x2).(x2-x1) / |x2-x1| . (x3-x2) / |x3-x2|^3
%        -(x3-x2) / |x4-x3| / |x3-x2|
%        +(x4-x3) / |x4-x3| / |x3-x2|
%        -(x4-x3).(x3-x2) / |x4-x3|^3 . (x3-x4)/|x3-x2|
%        -(x4-x3).(x3-x2) / |x4-x3| . (x3-x2)/|x3-x2|^3
%        -(x5-x4)/|x5-x4|/|x4-x3|
%        -(x5-x4).(x4-x3)/|x5-x4| . (x3-x4)/|x4-x3|^3
%      ]
%
% force on particle np has 3 contributions, from th_np, th_np+1, th_np+2
%
%
%
% force on particle Np-1 has 2 contributions, from th_Np-3 and th_Np-2
%
% force on particle Np has 1 contribution, from th_Np-2

for nt=1:Nt-1
    
    % elastic
    bR=norm(xp(:,1,nt)-xp(:,2,nt));
    f(:,1)=-k*(bR-b0)*(xp(:,1,nt)-xp(:,2,nt))/bR;
    for np=2:Np-1
        bL=norm(xp(:,np,nt)-xp(:,np-1,nt));
        fL=-k*(bL-b0)*(xp(:,np,nt)-xp(:,np-1,nt))/bL;
        bR=norm(xp(:,np,nt)-xp(:,np+1,nt));
        fR=-k*(bR-b0)*(xp(:,np,nt)-xp(:,np+1,nt))/bR;
        f(:,np)=fL+fR;  
                
        
    end
    bL=norm(xp(:,Np-1,nt)-xp(:,Np,nt));
    f(:,Np)=-k*(bL-b0)*(xp(:,Np,nt)-xp(:,Np-1,nt))/bL;
    
    % bending
    xc=xp(:,2,nt)-xp(:,1,nt);rc=abs(xc);
    xd=xp(:,3,nt)-xp(:,2,nt);rd=abs(rd);
    f(:,1)=f(:,1)+kb*(-xd/rd/rc + dot(xd,xc)/rd/rc^3 * xc);
    
    xb=xp(:,2,nt)-xp(:,1,nt);rb=abs(xb);
    xc=xp(:,3,nt)-xp(:,2,nt);rc=abs(xc);
    xd=xp(:,4,nt)-xp(:,3,nt);rd=abs(xd);
    f(:,2)=f(:,2)+kb*(-xb/rc/rb + xc/rc/rb                              ...
                      -dot(xc,xb)/rc/rb^3 * xb                          ...
                      +dot(xc,xb)/rb/rc^3 * xc                          ...
                      -xd/rd/rc + dot(xd,xc)/rd/rc^3 * xc);
                  
    for np=3:Np-2
        
        xa=xp(:,np-1,nt)-xp(:,np-2,nt);ra=abs(xa);
        xb=xp(:,np,nt)  -xp(:,np-1,nt);rb=abs(xb);
        xc=xp(:,np+1,nt)-xp(:,np,nt)  ;rc=abs(xc);
        xd=xp(:,np+2,nt)-xp(:,np+1,nt);rd=abs(xd);
        
        f(:,np)=f(:,np)+kb*( xa/ra/rb - dot(xb,xa)/rb^3/ra * xb         ...
                            -xb/rc/rb + xc/rc/rb                        ...
                            +dot(xc,xb)/rc^3/rb*xc - dot(xc,xb)/rc/rb^3 ...
                            -xd/rd/rc + dot(xd,xc)/rd/rc^3 * xc );
    end
    
    
    xa=xp(:,Np-2,nt)-xp(:,Np-3,nt);ra=abs(xa);
    xb=xp(:,Np-1,nt)-xp(:,Np-2,nt);rb=abs(xb);
    xc=xp(:,Np  ,nt)-xp(:,Np-1,nt);rc=abs(xc);
    
    f(:,Np-1)=f(:,Np-1)+kb*( xa/rb/ra - dot(xb,xa)/rb^3/ra * xb         ...
                            -xb/rc/rb + xc/rc/rb                        ...
                            +dot(xc,xb)/rc^3/rb * xc                    ...
                            -dot(xc,xb)/rc/rb^3 * xb );
    
    xa=xp(:,Np-1,nt)-xp(:,Np-2,nt);ra=abs(xa);
    xb=xp(:,Np  ,nt)-xp(:,Np-1,nt);rb=abs(xb);
    
    f(:,Np)=f(:,Np)+kb*( xa/rb/ra - dot(xb,xa)/rb^3/ra * xb);
                        
    % apply force to particles
    for np=1:Np
        f(:,np)=f(:,np)+sqrt(2*kB*T*R/dt)*randn(3,1);
        u=f(:,np)/R;
        xp(:,np,nt+1)=xp(:,np,nt)+u*dt;
    end    
end

figure(1);clf;hold on;
for np=1:Np
    plot3(squeeze(xp(1,np,:)),squeeze(xp(2,np,:)),squeeze(xp(3,np,:)));
end
plot3(squeeze(xp(1,:,1)), squeeze(xp(2,:,1)), squeeze(xp(3,:,1)), 'bo-');
plot3(squeeze(xp(1,:,Nt/2)),squeeze(xp(2,:,Nt/2)),squeeze(xp(3,:,Nt/2)),'ro-');
plot3(squeeze(xp(1,:,Nt)),squeeze(xp(2,:,Nt)),squeeze(xp(3,:,Nt)),'ko-');
% next steps... add brownian motion. then several particles, then several particles with bending
% potential. then string of particles.