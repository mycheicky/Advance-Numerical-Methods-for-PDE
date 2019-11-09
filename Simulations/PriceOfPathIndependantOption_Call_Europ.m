%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     PARAMETERS                      %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
clear;
clc;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     VARIABLES                       %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%Initialization parameters
T = 5;
r = 0.7;
sig = 0.3;
N = 41;
Nmc = 100;
At = T/N;
t = linspace(0,T,N+1);

S(1) = 10;
K = 10;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%                  Call Europ & Butterfly             %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% figure
% for n = 1:Nmc
%     for i = 1:N
%         S(i+1) = S(i) + S(1)*exp((r-(sig^2)/2)*At + sig*sqrt(T-t(i+1))*randn(1,1));
%     end
%     hold on
%     plot(t,S);
% end

GraphCall();

SurfCall();
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function [price] = PayOffCallO(So, T)
    gain = 0;
    Nmc = 10000;
    sig = 0.3;
    r = 0.07;
    t = linspace(0,T,Nmc+1);
    K = 10;
    for k = 1:Nmc
        ST = So*exp(r-(sig^2/2)*(T-t(k)) + sig*sqrt(T-t(k))*randn(1,1));
        gain = gain + max(ST - K, 0);
    end
    price = exp(-r*(T-t(k)))*gain/Nmc;
end

function [price] = PayOffCallButterfly(So, T)
    gain = 0;
    Nmc = 10000;
    sig = 0.3;
    r = 0.07;
    t = linspace(0,T,Nmc+1);
    K = 10;
    for k = 1:Nmc
        ST = So*exp(r-(sig^2/2)*(T-t(k)) + sig*sqrt(T-t(k))*randn(1,1));
        
        if ST > K && ST <= 2*K
            P = ST-K;
        elseif ST > 2*K && ST <= 3*K
            P = 3*K - ST;
        else
            P = 0;
        end
        
        gain = gain + P;
    end
    price = exp(-r*(T-t(k)))*gain/Nmc;
end

function[] = GraphCall()
    for i = 1:80
        So(i) = 0.5*(i-1);
        priceE(i) = PayOffCallO(So(i),5); 
        priceB(i) = PayOffCallButterfly(So(i),5);
    end
    
    figure
    plot(So,priceE);
    title('European Call Pay Off');
    xlabel('So');
    ylabel('Pay Off');
    
    figure
    plot(So,priceB);
    title('European Butterfly Pay Off');
    xlabel('So');
    ylabel('Pay Off');
end

function[] = SurfCall()
    for n = 1:21
        for i = 1:80
            St(i) =0.5*(i-1);
            time(n) = 0.025*(n-1);
            V_Eu_Call(n,i) = PayOffCallO(St(i),time(n));
            V_But_Call(n,i) = PayOffCallButterfly(St(i),time(n)); 
        end
    end
    
    %European Call
    figure
    surf(St,time,V_Eu_Call);
    title('European Call Surface');
    xlabel('St');
    ylabel('Time');
    zlabel('Call');
    
    %Butterfly
    figure
    surf(St,time,V_But_Call);
    title('European Butterfly Surface');
    xlabel('St');
    ylabel('Time');
    zlabel('Call');
    
end
