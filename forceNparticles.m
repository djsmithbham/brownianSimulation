% elastic force on 2 particles connected with a spring. Hydrodynamic drag
% but no hydrodynamic interaction

% SI
a=1e-6;
mu=1e-3;
k=2e-7; % not sure about this! spring constant
b0=3e-6; % equilibrium bond length
R=6*pi*mu*a;

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

% elastic potential is   1/2 k (b-b0)^2
%              = 1/2 k * 2 * (b-b0) * db/dxj
% but b = sqrt((xa-xb)^2+(ya-yb)^2+(za-zb)^2)
% so  db/dxj = 1/2 * 2 (xaj-xbj) / b
%            = (xaj-xbj) / b
% so force is - k (b-b0) (xaj-xbj) / b

% for a string of particles, have force from 'left' and force from 'right'
for nt=1:Nt-1
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
    for np=1:Np
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