function result=Heston_SZ4(K,delta,exact,scheme)
%delta: size of time step
%exact: exact value of the option
%scheme: 4:EAE,1:Euler,2:Euler_Central,3:BK_Euler
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

c0=-1/2*rho_xv*tau*delta;
c1=-delta1*rho_xv*psi*k*delta/tau;
c2=-delta2*rho_xv*psi*k*delta/tau;
c3=-1/2*delta1*delta+rho_xv/tau*(delta1*k*delta-1/2);
c4=-1/2*delta2*delta+rho_xv/tau*(delta2*k*delta+1/2);
c5=sqrt(1-rho_xv^2)*sqrt(delta);

d1=c1;
d2=c2;
d3=c3+1/2*(1-rho_xv^2)*delta1*delta;
d4=c4+1/2*(1-rho_xv^2)*delta2*delta;

e0=1/2*log(1-2*d4*k3^2)-d4*(k2+(d2/2/d4)^2)/(1-2*d4*k3^2);
e1=-d1-2*d4*k1*(k2+(d2/2/d4))/(1-2*d4*k3^2);
e2=-d3-d4*k1^2/(1-2*d4*k3^2);




N=time/delta;
Y=NaN(MM,1);

tic 
for l=1:MM
        for i=1:N
            if i==1
                v1=zeros(M,1);
                xx1=log(s0)*ones(M,1);
            end
            
            temp=randn(M,1)
            v2=k1*v1+k2+k3*temp;
            
         
            switch scheme
                case 4
                    c0_new=e0+e1*v1+e2*v1.^2;
                    xx2=xx1+c0_new+c1*v1+c2*v2+c3*v1.^2+c4*v2.^2+sqrt(delta1*v1.^2+delta2*v2.^2).*c5.*randn(M,1);
                case 1
                    xx2=xx1-1/2*v1.^2.*delta+v1.*sqrt(delta).*(rho_xv*temp+sqrt(1-rho_xv^2)*randn(M,1));
                case 2
                    xx2=xx1-1/2*(delta1*v1.^2+delta2*v2.^2)*delta+rho_xv*v1.*sqrt(delta).*randn(M,1)+...
                        sqrt(1-rho_xv^2)*sqrt(delta1*v1.^2+delta2*v2.^2).*sqrt(delta).*randn(M,1);
                case 3
                    c0_new=e0_3+e1_3*v1+e2_3*v1.^2;
                    xx2=xx1+c0_new+c1_3*v1+c2_3*v2+c3_3*v1.^2+c4_3*v2.^2+v1.*c5_3.*randn(M,1);
            end
        
            v1=v2;
            xx1=xx2; 
        end
     
    Y(l)=mean(max(exp(xx1)-K,0));
end
toc
save('Y.mat','y_temp');

Y=Y-exact;
save('Y.mat','y');
result=[mean(Y);std(Y)*2.326]; % 99% CI

end