function result=Heston_Euler(K,delta)
%Euler
if (K==100)
    exact = 27.9;
elseif (K==140)
    exact = 14.23;
elseif (K==60)
    exact = 50.34;
else
    fprintf('error K\n');
end

k=0.1;
tau=0.3;
v0=0;
s0=100;
theta=0;
rho_xv=-0.6;
r=0;
psi=0;
time=5; %option time length
M=1000000; %number of simulation paths
MM=1; %number of trials

k1=exp(-k*delta);
k2=psi*(1-exp(-k*delta));
k3=tau*sqrt((1-exp(-2*k*delta))/2/k);

N=time/delta;
Y=NaN(MM,1);
SD=Y;

tic
for l=1:MM
    v1=zeros(M,1);
    xx1=log(s0)*ones(M,1);
    for i=1:N
        temp=randn(M,1);
        v2=k1*v1+k2+k3*temp;
        xx2=xx1-1/2*v1.^2.*delta+v1.*sqrt(delta).*(rho_xv*temp+sqrt(1-rho_xv^2)*randn(M,1));
        
        v1=v2;
        xx1=xx2;
    end
    payoff = max(exp(xx1)-K,0);
    Diff = payoff-exact;
    Y(l)=mean(Diff);
    SD(l)=std(Diff)/sqrt(M)*2.326;
end
toc
% save('Y_temp.mat','Y');
% 
% YY=Y-exact;
% save('Y.mat','YY');
% result=[mean(YY);std(YY)*2.326/sqrt(MM)]; % 99% CI
result=[Y SD];
end