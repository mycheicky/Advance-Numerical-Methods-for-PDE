%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     PARAMETERS                      %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
clear;
clc;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                     VARIABLES                       %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%Initialization parameters
T = 1;
u = 0.12;
p = 2;
r = 0.07;
Sig = [0.3;0.4];
lam = 0;
N = 100;
Nmc = 100000;
At = T/N;
t = linspace(0,T,N+1);

count = 0;

X(1) = 100;
M(1) = 500;
P(1) = 100;
Oo = 0:0.1:1;
FkExp = 0;
j = 1;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%         Underlying Asset Monte Carlo Simulation     %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

figure
for n = 1:Nmc
    for i = 1: N
        X(i+1) = X(i) - p*X(i)*At + Sig(2)*At*randn;
        M(i+1) = M(i) + M(i)*r*At*(1-Oo(j)) + M(i)*((u+lam*X(i))*At + Sig(1)*randn*sqrt(At))*Oo(j);
        P(i+1) = P(i) + ((u+lam*X(i))*At + Sig(1)*randn*sqrt(At))*P(i);

        if M(i) < 800
            count = count + 1;
        end
    end

    MT(n) = M(N+1);

    hold on
    plot(t,M);
    
    if j < length(Oo)
        j = j + 1;
    end
end
title('Graph of Mt');
ylabel('Mt');
xlabel('t');
% 
% PMt = count/Nmc;
% Exp(j) = sum(log(MT))/Nmc;

% figure
% plot(t,X);
% title('Graph of X');
% ylabel('Xt');
% xlabel('t');
% 
% figure
% plot(t,P);
% title('Graph of P');
% ylabel('Pt');
% xlabel('t');
% 
% figure
% plot(Oo,Exp);
% title('Graph of Exp');
% ylabel('Exp');
% xlabel('t');
% 
% disp('Probability =');
% disp(PMt);
% 
% disp('Expectation Value =');
% disp(max(Exp));
% 
% disp('Theta Value =');
% disp(MaxTheta);

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

