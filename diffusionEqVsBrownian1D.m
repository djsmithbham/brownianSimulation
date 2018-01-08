clear all
close all

Nx=1000;
Nt=500;

x=linspace(-0.5,0.5,Nx);
t=linspace(0,0.05,Nt);

sig=sqrt(1/1000/2);

[X T]=ndgrid(x,t);

dx=x(2)-x(1);
dt=t(2)-t(1);

a=dt/dx^2;

U=zeros(Nx,Nt);

U(:,1)=exp(-x.^2/2/sig^2);
U(:,1)=U(:,1)/sum(U(:,1))/dx;

figure(1);clf;
subplot(2,2,1);plot(x,U(:,1));

e = ones(Nx,1);
A = spdiags([-a/2*e (1+a)*e -a/2*e], -1:1, Nx, Nx);
B = spdiags([ a/2*e (1-a)*e  a/2*e], -1:1, Nx, Nx);
A(1,1)  =1; A(1,2)    =-1;
A(Nx,Nx)=1; A(Nx,Nx-1)=-1;
B(1,1)  =0; B(1,2)    =0;
B(Nx,Nx)=0; B(Nx,Nx-1)=0;

for k=1:Nt-1
    U(:,k+1)=A\(B*U(:,k));
end

subplot(2,2,2);pcolor(X,T,U);shading flat;colorbar;

subplot(2,2,3);plot(t,sum(U,1));xlabel('t');ylabel('\int_{-L/2}^{L/2} U(x,t) dx','interpreter','tex');

%%

Nsamp=5000;

xp=zeros(Nt,Nsamp);
xp(1,:)=sig*randn([1,Nsamp]);

mean(xp(1,:))
std(xp(1,:))

figure(1);
subplot(2,2,4);h=histogram(xp(1,:),'Normalization','pdf');xlim([-0.5 0.5]);

for s=1:Nsamp
    for nt=1:Nt-1
        xp(nt+1,s)=xp(nt,s)+sqrt(2*dt)*randn;
    end
end

nt=floor(Nt/8);

figure(2);clf;

subplot(1,2,1);plot(x,U(:,nt));ylim([0 5]);

subplot(1,2,2);h=histogram(xp(nt,:),'Normalization','pdf');xlim([-0.5,0.5]);ylim([0 5]);


nt=floor(Nt/4);

figure(3);clf;

subplot(1,2,1);plot(x,U(:,nt));ylim([0 5]);

subplot(1,2,2);h=histogram(xp(nt,:),'Normalization','pdf');xlim([-0.5,0.5]);ylim([0 5]);