clear;
clc;
T=0.5;
r=0.1;
sig=0.5;
t=0;
S0=10;
K=10;
Nmc=10;
Tt=1.73;
theta=Tt;
sigma=sig;
Graph_Option_BS(K,T,r,sig)
for i=1:200
    
    S0(i)=0.1*(i-1);
    IST(i)=CallIST(Nmc,r,sigma,T,K,S0(i),theta);
end
hold on
plot(S0,IST,'r')
legend('Monte-Carlo','Importance Simpling Technic')
function [f]=CallIST(Nmc,r,sigma,T,K,S0,theta)
sum3=0
for i=1:Nmc
        WT=sqrt(T)*randn;
        g=rand;
        ST2=S0* exp(((r-((sigma^2)/2))*T) + (sigma * (WT + theta*T)));
        sum3=sum3+ Function_Payoff_Call(ST2,K)* exp(-theta*WT - (theta*theta*T)/2);  
end
f=sum3*exp(-r*T)/Nmc
end


function []=Graph_Option_BS(K,T,r,sig)
for j=1:41
    St(j)= 0.5*(j-1);
    price(j)=BS_theory1(St(j),r,T,K,sig);
end
% hold on
% plot(St,price,'r')
% title('Pay-Off at t=0 ')
% xlabel('Underlying Asset')
% ylabel('Price')
% grid on
end
function[f]=Function_Payoff_Call(S,K)

f=max(S-K,0);
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