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
N = 41;
Nmc = 1000;
At = T/N;
t = linspace(0,T,N+1);

So = 10;
K = 10;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%                  Call Europ & Butterfly             %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%Values for the confidence interval
C = Callo(So,K,r,T,sig,N,Nmc);
disp('Expectation =');
disp(C(1));
disp('Variance =');
disp(C(2));
disp('Lower value =');
disp(C(3));
disp('Upper value =');
disp(C(4));

%Probability of the confidence interval
P = ConfInterProba(So,K,r,T,sig,N,Nmc);
disp('Probability of the confident interval =');
disp(P);



%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [WT] = Q(So,K,r,T,sig,g,N)
    WT = exp(-r*T)*max(So * exp((r-sig^2/2)*T + sig*g*sqrt(T)) - K,0);
end

function [C] = Callo(So,K,r,T,sig,N,Nmc)
    for i = 1:Nmc
        %g = randn(1,1);
        E(i) = Q(So,K,r,T,sig,randn(1,1),N);
        %exp(-r*T)*max(So * exp((r-(sig*sig)/2)*T + sig*randn(1,1)*sqrt(T)) - K,0);
        E2(i) = Q(So,K,r,T,sig,randn(1,1),N)^2;
        %( exp(-r*T)*max(So * exp((r-sig^2/2)*T + sig*randn(1,1)*sqrt(T)) - K,0) )^2;
    end
    Exp = sum(E)/Nmc;
    Var = sum(E2)/Nmc - Exp^2;
    %The Confidence Interval
    lower =  Exp - 1.96*sqrt(Var)/sqrt(Nmc);
    upper =  Exp + 1.96*sqrt(Var)/sqrt(Nmc);
    C = [Exp,Var,lower,upper];
end

function [P] = ConfInterProba(So,K,r,T,sig,N,Nmc)
    Proba = 0;
    for i = 1:Nmc
        C = Callo(So,K,r,T,sig,N,Nmc);
        if C(3) < C(1) && C(1) < C(4)
            Proba = Proba + 1;
        end
    end
    P = Proba/Nmc;
end
