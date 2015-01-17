function result=SZ_4(K,delta,exact)
%EAE scheme
k=0.4;
tau=0.4;
v0=0.2;
s0=100;
%theta=0.2;
rho_xv=-0.9;
r=0.04;
psi=0.2;
time=10; %option time length
M=10000; %number of simulation paths
MM=1000; %number of trials

delta1=1/2;
delta2=1/2;

k1=exp(-k*delta);
k2=psi*(1-exp(-k*delta));
k3=tau*sqrt((1-exp(-2*k*delta))/2/k);

%c0=-1/2*rho_xv*tau*delta;
c1=-delta1*rho_xv*psi*k*delta/tau;
c2=-delta2*rho_xv*psi*k*delta/tau;
c3=-1/2*delta1*delta+rho_xv/tau*(delta1*k*delta-1/2);
c4=-1/2*delta2*delta+rho_xv/tau*(delta2*k*delta+1/2);
c5=sqrt(1-rho_xv^2)*sqrt(delta);

d1=c1;
d2=c2;
d3=c3+1/2*(1-rho_xv^2)*delta1*delta;
d4=c4+1/2*(1-rho_xv^2)*delta2*delta;

e0=1/2*log(1-2*d4*k3^2)-d4*(k2+(d2/2/d4))^2/(1-2*d4*k3^2);
e1=-d1-2*d4*k1*(k2+(d2/2/d4))/(1-2*d4*k3^2);
e2=-d3-d4*k1^2/(1-2*d4*k3^2);


N=time/delta;
Y=NaN(MM,1);

tic
for l=1:MM
    v1=ones(M,1)*v0;
    xx1=log(s0)*ones(M,1);
    for i=1:N
        v2=k1*v1+k2+k3*randn(M,1);
        c0_new=e0+e1*v1+e2*v1.^2;
        xx2=xx1+c0_new+c1*v1+c2*v2+c3*v1.^2+c4*v2.^2+sqrt(delta1*v1.^2+delta2*v2.^2).*c5.*randn(M,1);
        
        v1=v2;
        xx1=xx2;
    end
    Y(l)=mean(max(exp(xx1)-K,0));
end
toc
save('Y_temp.mat','Y');

YY=Y-exact;
save('Y.mat','YY');
result=[mean(YY);std(YY)*2.326]; % 99% CI
end