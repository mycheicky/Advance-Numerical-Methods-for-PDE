%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     PARAMETERS                      %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
clear;
clc;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     VARIABLES                       %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%Initialization parameters
T = 0.5;
r = 0.1;
sig = 0.5;
N = 100;
Nmc = 10000;
At = T/N;
t = linspace(0,T,N+1);

So = 10;
K = 10;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% UnderP(10000000);
% UnderQ(1000,5);

theta = 1.782;
C2 = Call3(So,K,r,T,theta,sig,N,Nmc);
C0 = Callo(So,K,r,T,sig,N,Nmc);
disp('Expectation Call0 =');
disp(C0(1));
disp('Expectation Call3 =');
disp(C2(1));
disp('Call0 =');
disp(C0(2));
disp('Call3 =');
disp(C2(2));

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                      FUNCTIONS                      %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [] = UnderP(Nmc)
    count = 0;
    for i = 1:Nmc
        g = randn(1,1);
        if g >=5
            count = count + 1;
        end
    end
    P = count/Nmc;
    disp('Under P we got IP =');
    disp(P);
end

function [] = UnderQ(Nmc,miu)
    count = 0;
    for i = 1:Nmc
        g = randn;
        if g + miu >=5
            count = count + exp(-g*miu - (miu^2)/2);
        end
    end
    P = count/Nmc;
    
    disp('Under p -> Q we got IP =');
    disp(P);
end

function [WT] = Qw(So,K,r,T,sig,g,N)
    WT = exp(-r*T)*max(So * exp((r-sig^2/2)*T + sig*g*sqrt(T)) - K,0);
end

function [C] = Callo(So,K,r,T,sig,N,Nmc)
    for i = 1:Nmc
        g = randn(1,1);
        E(i) = Qw(So,K,r,T,sig,g,N);
        E2(i) = Qw(So,K,r,T,sig,g,N)^2;
    end
    
    Exp = sum(E)/Nmc;
    Var = sum(E2)/Nmc - Exp^2;
    
    %The Confidence Interval
    lower =  Exp - 1.96*sqrt(Var)/sqrt(Nmc);
    upper =  Exp + 1.96*sqrt(Var)/sqrt(Nmc);
    
    C = [Exp,Var,lower,upper];
end

function [WT] = Q(So,K,r,T,theta,sig,g,N)
    wt = g + theta*T;
    WT = exp(-r*T)*max(So * exp((r-sig^2/2)*T + sig*wt*sqrt(T)) - K,0)*exp(-theta*wt - (theta^2)*T/2);
end

function [C] = Call3(So,K,r,T,theta,sig,N,Nmc)
    for i = 1:Nmc
        g = randn(1,1);
        E(i) = Q(So,K,r,T,theta,sig,g,N);
        E2(i) = Q(So,K,r,T,theta,sig,g,N)^2;
    end
    
    Exp = sum(E)/Nmc;
    Va = sum(E2)/Nmc - Exp^2;
    
    %The Confidence Interval
    lower =  Exp - 1.96*sqrt(Va)/sqrt(Nmc);
    upper =  Exp + 1.96*sqrt(Va)/sqrt(Nmc);
    
    C = [Exp,Va,lower,upper];
end
