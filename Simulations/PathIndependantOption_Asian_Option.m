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
r = 0.4;
sig = 0.5;
N = 100;
Nmc = 1000;
% At = 0.5;
% t = linspace(0,T,N+1);

So = 10;
K = 10;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%                  Call Europ & Butterfly             %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%
figure
GraphAsian(0,10,r,sig,T,N,Nmc,K,50);

% figure
% SurfAsian(T/3,r,sig,T,N,Nmc,K);



%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%%% Graph %%%

%
function[payoff] = AsianOpt_PO_1(t,S0,r,sig,T,N,K)
    At = T/N;
    A = 0;
    S(1) = S0;
    for i = 1:N
        S(i+1) = S(i)*exp( (r-(sig^2)/2)*At + sig*sqrt(At)*randn(1,1) );
        A = A + S(i+1)*At;
    end
    payoff = max(A/T - K,0);
end

%
function[payoff] = EuroCallOpt_PO_1(t,S0,r,sig,T,N,K)
    At = T/N;
    A = 0;
    S(1) = S0;
    for i = 1:N
        S(i+1) = S(i)*exp( (r-(sig^2)/2)*At + sig*sqrt(At)*randn(1,1) );
    end
    payoff = max(S(i+1) - K,0);
end

%
function[price] = PriceAsianOpt_S0_1(t,St,At,r,sig,T,N,Nmc,K)
    P = 0;
    for i = 1:Nmc
        P = P + AsianOpt_PO_1(t,St,r,sig,T,N,K);
    end
    price = exp(-r*T) * P/Nmc;
end

%
function[price] = PriceEuroCallOpt_S0_1(t,St,At,r,sig,T,N,Nmc,K)
    P = 0;
    for i = 1:Nmc
        P = P + EuroCallOpt_PO_1(t,St,r,sig,T,N,K);
    end
    price = exp(-r*T) * P/Nmc;
end

%
function[] = GraphAsian(t,S0,r,sig,T,N,Nmc,K,M)
    At = 2*K/M;
    for i = 1:M+1
        So(i) = At*(i-1);
        asian(i) = PriceAsianOpt_S0_1(0,So(i),At,r,sig,T,N,Nmc,K);
        europ(i) = PriceEuroCallOpt_S0_1(t,So(i),At,r,sig,T,N,Nmc,K);
        PO(i) = max(So(i) - K,0);
    end
    plot(So,europ,So,asian,So,PO);
    xlabel('Time (t)');
    ylabel('Asset (St)');
    legend('European Option','Asian Option','European Option Pay OffS');
    title('Call European Option versus Call Asian Option'),
    
end

%%% Surface %%%

%
function[payoff] = AsianOpt_PO_2(t,S0,At,r,sig,T,N,K)
    A = At*t;
    S(1) = S0;
    for i = 1:N
        S(i+1) = S(i)*exp( (r-(sig^2)/2)*At + sig*sqrt(At)*randn(1,1) );
        A = A + S(i+1)*At;
    end
    payoff = max(A/T-K,0);
end

%
function[price] = PriceAsianOpt_S0(t,St,At,r,sig,T,N,Nmc,K)
    P = 0;
    for i = 1:Nmc
        P = P + AsianOpt_PO_2(t,St,At,r,sig,T,N,K);
    end
    price = exp( -r*(T-t) ) * P/Nmc;
end

%
function[] = SurfAsian(t,r,sig,T,N,Nmc,K)
    for j = 1:41
        for n = 1:41
            St(j) = 0.5*(j-1);
            At(n) = 0.5*(n-1);
            price(n,j) = PriceAsianOpt_S0(t,St(j),At(n),r,sig,T,N,Nmc,K);
        end
    end
    surf(St,At,price);
end