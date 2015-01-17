function result=Heston_3(K,delta,exact)
%BK_Euler scheme
k=0.1;
tau=0.3;
v0=0;
s0=100;
theta=0;
rho_xv=-0.6;
r=0;
psi=0;
time=5; %option time length
M=10000; %number of simulation paths
MM=1000; %number of trials


k1=exp(-k*delta);
k2=psi*(1-exp(-k*delta));
k3=tau*sqrt((1-exp(-2*k*delta))/2/k);

%c0=-1/2*rho_xv*tau*delta;
c1=-delta*rho_xv*psi*k/tau;
c2=0;
c3=-1/2*delta+rho_xv/tau*(k*delta-1/2);
c4=1/2*rho_xv/tau;
c5=sqrt(1-rho_xv^2)*sqrt(delta);

d1=c1;
d2=c2;
d3=c3+1/2*(1-rho_xv^2)*delta;
d4=c4;

e0=1/2*log(1-2*d4*k3^2)-d4*(k2+(d2/2/d4))^2/(1-2*d4*k3^2);
e1=-d1-2*d4*k1*(k2+(d2/2/d4))/(1-2*d4*k3^2);
e2=d3-d4*k1^2/(1-2*d4*k3^2);

N=time/delta;
Y=NaN(MM,1);

tic
for l=1:MM
    v1=zeros(M,1);
    xx1=log(s0)*ones(M,1);
    for i=1:N
        v2=k1*v1+k2+k3*randn(M,1);
        c0_new=e0+e1*v1+e2*v1.^2;
        xx2=xx1+c0_new+c1*v1+c2*v2+c3*v1.^2+c4*v2.^2+v1*c5.*randn(M,1);
        
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