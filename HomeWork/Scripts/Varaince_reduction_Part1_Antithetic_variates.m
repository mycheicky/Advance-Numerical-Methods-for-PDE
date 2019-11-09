
clear;
clc;
close all;
T=0.5;
r=0.1;
sig=0.5;
t=0;
S0=10;
K=10;
Nmc=25;



% Graph_Option_BS(K,T,r,sig)

for i=1:200
    st(i)=0.1*(i-1);
    method1(i)=exp_val(Nmc,st(i),r,sig,T,t,K);
    var_red(i)=var_reduction(Nmc,st(i),r,sig,T,t,K);
end
hold on
plot(st,method1,'b')
% hold on
% plot(st,var_red,'r')
% grid on
%title('Monte Carlo expectation vs')
title('Method 2 : Control Variates')
legend('Monte-Carlo','Antithetic Variaties')

function [f]=exp_val(Nmc,S0,r,sig,T,t,K)
for j=1:Nmc
O(j)=exp(-r*T)*max(S0*exp((r-sig*sig/2)*(T-t)+sig*sqrt(T-t)*randn)-K,0);
end
exp_O=sum(O)/Nmc;
f=exp_O;
end

function [f]=var_reduction(Nmc,st,r,sig,T,t,K)
for j=1:Nmc
    g=randn;
    O1(j)=exp(-r*T)*max(st*exp((r-sig*sig/2)*(T-t)+sig*sqrt(T-t)*g)-K,0);
    O2(j)=exp(-r*T)*max(st*exp((r-sig*sig/2)*(T-t)+sig*sqrt(T-t)*(-g))-K,0);
end
Z=sum((O1+O2)/2)/Nmc;
f=Z;
end

function []=Graph_Option_BS(K,T,r,sig)
for j=1:41
    St(j)= 0.5*(j-1);
    price(j)=BS_theory1(St(j),r,T,K,sig);
end
hold on
plot(St,price,'k')
%title('pay-off at initial time ')
xlabel('Underlying Asset')
ylabel('Price')
grid on
end

function[f]=BS_theory1(S,r,T,K,sigma)
t=0;
f=S*N(d1(S,r,T,K,sigma))-K*exp(-r*(T-t))*N(d2(S,r,T,K,sigma));
end 

function[f]=d1(S,r,T,K,sigma)
t=0;
f=(log(S/K)+(r+sigma^2/2)*(T-t))/(sigma*sqrt(T-t));
end

function[f]=d2(S,r,T,K,sigma)
t=0;
f=(log(S/K)+(r-sigma^2/2)*(T-t))/(sigma*sqrt(T-t));
end

function[f]=N(x)
f=1/2*(1+erf(x/sqrt(2)));
end