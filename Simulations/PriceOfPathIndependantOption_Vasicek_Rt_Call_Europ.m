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

AllCurve();

% SurfCall();

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

function [price] = PayOffCallVasicek(So, T)
    gain = 0;
    ST(1) = So;
    Nmc = 10000;
    N = 100;
    At = T/N;
    sig = 0.3;
    r = Vasicek();
    t = (0:N)*At;
    K = 10;
    for i = 1:Nmc
        for k = 1:N
            ST = So*exp( r(k)-(sig^2/2)*At + sig*sqrt(At)*randn );
            gain = gain + exp(-r(k)*(T-t(k)))*max(ST - K, 0)/N;
        end
    end
    price = gain/Nmc;
end

function[Rt] = Vasicek()
    Nmc = 10000;
    r(1) = 0.01;
    gamma = 0.2;
    omega = 0.02;
    miu = 0.016;
    T = 5;
    N = 100;
    At = T/N;
    
    for i = 1:N
        r(i+1) = r(i)*(1-gamma*At) + miu*At + omega*randn*sqrt(At);
        
%         r(i+1)=r(i)*(1-gamma*delta_t)+nu*delta_t+omega*g*sqrt(delta_t);
%         r(i+1) = r(i) + miu/gamma + (r(i)-miu/gamma)*exp(-gamma*At) + omega*exp(-gamma*At)*randn*sqrt(At);
    end
    
    Rt = r;
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

function[] = AllCurve()
    for i = 1:80
        So(i) = 0.5*(i-1);
        price(i) = PayOffCallO(So(i),5); 
        priceV(i) = PayOffCallVasicek(So(i),5);
    end
    
    figure
    plot(So,price);
    title('European Call Pay Off');
    xlabel('Asset (S)');
    ylabel('Pay Off');
    
    figure
    plot(So,priceV);
    title('European Call Pay Off using Vasicek');
    xlabel('Asset (S)');
    ylabel('Pay Off');
    
    figure
    plot(So,price,So,priceV);
    title('Fixed Rate vs Vasicek Rate');
    xlabel('Asset (S)');
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
    
    plot(time,V_But_Call);
    
end
