clear all

Nr=1000;
Nt=500;

r=linspace(0,1,Nr);
t=linspace(0,0.02,Nt);

sig=sqrt(1/1000/2);

[R T]=ndgrid(r,t);

dr=r(2)-r(1);
dt=t(2)-t(1);

Nsamp=5000;

a=dt/dr^2;
b=dt/dr;

U=zeros(Nr,Nt);

U(:,1)=(20*r<1)*20; % step function
U(:,1)=filter(ones(3,1)/3,1,U(:,1));  % moving average filter
U(:,1)=U(:,1)/(sum(U(:,1).*r.^2')*4*pi*dr);

figure(1);clf;
subplot(2,2,1);plot(r,U(:,1).*r'.^2*4*pi);title(' 4 pi concentration * r^2, i.e. initial pdf');

e = ones(Nr,1);
A = spdiags([-a/2*e+b/2./[r(2:end) 1]' , (1+a)*e , -a/2*e-b/2./[1 r(1:end-1)]'], -1:1, Nr, Nr);
B = spdiags([ a/2*e-b/2./[r(2:end) 1]' , (1-a)*e ,  a/2*e+b/2./[1 r(1:end-1)]'], -1:1, Nr, Nr);
A(1,1)  =1; A(1,2)    =-1;
A(Nr,Nr)=1; A(Nr,Nr-1)=-1;
B(1,1)  =0; B(1,2)    =0;
B(Nr,Nr)=0; B(Nr,Nr-1)=0;

for k=1:Nt-1
    U(:,k+1)=A\(B*U(:,k));
end

subplot(2,2,2);pcolor(R,T,U);shading flat;colorbar;title('concentration');

subplot(2,2,3);plot(t,sum(U.*R.^2,1)*dr*4*pi);xlabel('t');ylabel('4\pi\int_{0}^{L} U(r,t) r^2 dr','interpreter','tex');title('check conservation for diffusion eq');

%%  

rp=zeros(Nt,Nsamp);
thp=zeros(Nt,Nsamp);

rp(1,:) =(rand([1,Nsamp])).^(1/3)/20; % produce uniform concentration in 
                                   % [0,0.05] (x^(1/3) is inverse of CDF
                                   % of p(x)=x^2)
                                   
thp(1,:)  = acos(1-2*rand([1,Nsamp]));    %  pi*rand([1,Nsamp]);
phip(1,:) =     2*pi*rand([1,Nsamp]);


figure(1);
subplot(2,2,4);h=histogram(rp(1,:),'Normalization','pdf');xlim([0 1]);title('simulation initial pdf');

xp(1,:)=rp(1,:).*sin(thp(1,:)).*cos(phip(1,:));
yp(1,:)=rp(1,:).*sin(thp(1,:)).*sin(phip(1,:));
zp(1,:)=rp(1,:).*cos(thp(1,:));


for s=1:Nsamp
    for nt=1:Nt-1
        xp(nt+1,s)=xp(nt,s)+sqrt(2*dt)*randn;
        yp(nt+1,s)=yp(nt,s)+sqrt(2*dt)*randn;
        zp(nt+1,s)=zp(nt,s)+sqrt(2*dt)*randn;
    end
end

rp=sqrt(xp.^2+yp.^2+zp.^2);

%% this section is not finished

figure(4);clf;hold on;
plot3(xp(:),yp(:),zp(:),'.');

%[thp,rp]=cart2pol(xp,yp);

nt=floor(Nt/8);

figure(2);clf;

subplot(1,3,1);plot(r,4*pi*U(:,nt).*r.^2');title('4 pi * concentration * r^2, i.e. pdf');

subplot(1,3,2);h=histogram(rp(nt,:),'Normalization','pdf');xlim([0,1]);title('Brownian simulation pdf');

subplot(1,3,3);plot3(xp(nt,:),yp(nt,:),zp(nt,:),'.');title('Brownian simulation (x,y,z)');


nt=floor(Nt/4);

figure(3);clf;

subplot(1,3,1);plot(r,4*pi*U(:,nt).*r.^2');title('4 pi * concentration * r^2, i.e. pdf');

subplot(1,3,2);h=histogram(rp(nt,:),'Normalization','pdf');xlim([0,1]);title('Brownian simulation pdf');

subplot(1,3,3);plot3(xp(nt,:),yp(nt,:),zp(nt,:),'.');title('Brownian simulation (x,y,z)');
