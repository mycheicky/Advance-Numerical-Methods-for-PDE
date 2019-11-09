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
Nmc = 25;
At = T/N;
theta = 1.73;

So = 10;
K = 10;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%                  Call Europ & Butterfly             %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%Graph
for i=1:200
    S(i)=0.1*(i-1);
    
    %Antithetic Variates
    C = VarianceReductionAV(S(i),K,r,T,sig,N,Nmc);
    V(i) = C(1); %Expectation
    AV(i) = C(2); %Variance Reduction
    
    %Control Variates
    CV(i) = VarianceReductionCV(S(i),K,r,T,sig,N,Nmc);
    
    %Importance Sampling Technique
    IST(i) = VarianceReductionIST(S(i),K,r,T,theta,sig,N,Nmc);
end

% Graph 1 : Monte Carlo Expectation
figure
plot(S,V,'b');
xlabel('Option Price (Vt)');
ylabel('Asset (St)');
title('Monte Carlo Expectation')
legend('Monte-Carlo Expectation')

% Graph 2 : Monte Carlo Expectation VS Antithetic Variates
figure
hold on
plot(S,V,'b');
hold on
plot(S,AV,'r');
grid on
xlabel('Option Price (Vt)');
ylabel('Asset (St)');
title('Method 1 : Antithetic Variates');
legend('Monte-Carlo','Antithetic Variaties');

% Graph 3 : Monte Carlo Expectation VS Control Variates
figure
hold on
plot(S,V,'b');
hold on
plot(S,CV,'r');
grid on
xlabel('Option Price (Vt)');
ylabel('Asset (St)');
title('Method 2 : Control Variates');
legend('Monte-Carlo','Control Variaties');

% Graph 4 : Monte Carlo Expectation VS Importance Sampling Technique
figure
hold on
plot(S,V,'b');
hold on
plot(S,IST,'r');
grid on
xlabel('Option Price (Vt)');
ylabel('Asset (St)');
title('Method 3 : Importance Sampling Technique');
legend('Monte-Carlo','Importance Sampling Technique');


%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [WT] = Qw(So,K,r,T,sig,g,N)
    WT = exp(-r*T)*max(So * exp((r-(sig^2)/2)*T + sig*g*sqrt(T)) - K,0);
end

function [C] = VarianceReductionAV(So,K,r,T,sig,N,Nmc) %Antithetic Variates
    for i = 1:Nmc
        g = randn(1,1);
        E(i) = Qw(So,K,r,T,sig,g,N); %Monte Carlo Expectation Value
        E2(i) = Qw(So,K,r,T,sig,-g,N);
    end
    Exp = sum(E)/Nmc;
    Var = sum( (E+E2)/2 )/Nmc;
    
    %The Confidence Interval
    lower =  Exp - 1.96*sqrt(Var)/sqrt(Nmc);
    upper =  Exp + 1.96*sqrt(Var)/sqrt(Nmc);
    
    C = [Exp,Var,lower,upper];
end

function [f] = VarianceReductionCV(So,K,r,T,sig,N,Nmc) %Control Variates => value of b
    for i = 1:Nmc
        St(i) = So * exp( ((r-((sig^2)/2))*T) + (sig * sqrt(T) * randn(1,1)) );
        Cov(i) = (exp(-r*T)*max(St(i)-K,0) - Qw(So,K,r,T,sig,randn(1,1),N)) * (St(i) - So*exp(r*T));
        Var(i) = (St(i) - So*exp(r*T))^2;
    end
    
    b = sum(Cov)/sum(Var);
    
    for i = 1:Nmc
        CV(i) = (exp(-r*T)*max(St(i)-K,0) - b * (St(i) - So*exp(r*T)));
    end
    
    f = sum(CV)/Nmc;
end

function [f] = W(So,K,r,T,theta,sig,g,N)
    wt = sqrt(T)*g;
    f = max(So * exp((r-(sig^2)/2)*T + sig*(wt + theta*T)) - K,0) * exp(-theta*wt - (theta^2)*T/2);
end

function [C] = VarianceReductionIST(So,K,r,T,theta,sig,N,Nmc)
    for i = 1:Nmc
        E(i) = W(So,K,r,T,theta,sig,randn(1,1),N);
    end
    
    Var = (sum(E) * exp(-r*T))/Nmc;
    
    C = Var;
end

function [P] = ConfInterProba(So,K,r,T,sig,N,Nmc)
    Proba = 0;
    for i = 1:Nmc
        C = VarianceReductionAV(So,K,r,T,sig,N,Nmc);
        if C(3) < C(1) && C(1) < C(4)
            Proba = Proba + 1;
        end
    end
    P = Proba/Nmc;
end
