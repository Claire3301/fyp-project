function result=Heston_2(K,delta,exact)
%Euler-central
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


delta1=1/2;
delta2=1/2;

k1=exp(-k*delta);
k2=psi*(1-exp(-k*delta));
k3=tau*sqrt((1-exp(-2*k*delta))/2/k);

N=time/delta;
Y=NaN(MM,1);

tic
for l=1:MM
    v1=zeros(M,1);
    xx1=log(s0)*ones(M,1);
    for i=1:N
        temp=randn(M,1);
        v2=k1*v1+k2+k3*temp;
        xx2=xx1-1/2*(delta1*v1.^2+delta2*v2.^2)*delta+rho_xv*v1.*sqrt(delta).*temp+...
        sqrt(1-rho_xv^2)*sqrt(delta1*v1.^2+delta2*v2.^2).*sqrt(delta).*randn(M,1);
        
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