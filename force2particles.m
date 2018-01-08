% elastic force on 2 particles connected with a spring. Hydrodynamic drag
% but no hydrodynamic interaction

% SI
a=1e-6;
mu=1e-3;
k=1e-7; % not sure about this! spring constant
b0=3e-6; % equilibrium bond length
R=6*pi*mu*a;

Np=2;
Nt=100;
tMax=1;

t=linspace(0,tMax,Nt);
dt=t(2)-t(1);

xp=zeros(3,Np,Nt);

xp(1,2,1)=5*a;
f=zeros(3,Np);

% elastic potential is   1/2 k (b-b0)^2
%              = 1/2 k * 2 * (b-b0) * db/dxj
% but b = sqrt((xa-xb)^2+(ya-yb)^2+(za-zb)^2)
% so  db/dxj = 1/2 * 2 (xaj-xbj) / b
%            = (xaj-xbj) / b
% so force is - k (b-b0) (xaj-xbj) / b

for nt=1:Nt-1
    b=norm(xp(:,1,nt)-xp(:,2,nt));
    f(:,1)=-k*(b-b0)*(xp(:,1,nt)-xp(:,2,nt))/b;
    f(:,2)=-k*(b-b0)*(xp(:,2,nt)-xp(:,1,nt))/b;
    for np=1:Np
        u=f(:,np)/R;
        xp(:,np,nt+1)=xp(:,np,nt)+u*dt;
    end
end

figure(1);clf;
plot(squeeze(xp(1,1,:)),'b');
hold on;
plot(squeeze(xp(1,2,:)),'r');


% next steps... add brownian motion. then several particles, then several particles with bending
% potential. then string of particles.