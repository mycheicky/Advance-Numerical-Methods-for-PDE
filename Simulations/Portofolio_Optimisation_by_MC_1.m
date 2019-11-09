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

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

MT(1) = M(1);
U(1) = log(MT(1));
figure
for j = 1:length(Oo) 
    for n = 1:Nmc
        for i = 1: N
            X(i+1) = X(i) - p*X(i)*At + Sig(2)*sqrt(At)*randn;
            M(i+1) = M(i) + M(i)*r*At*(1-Oo(j)) + M(i)*( (u+lam*X(i))*At + Sig(1)*randn*sqrt(At) )*Oo(j);
            P(i+1) = P(i) + ((u+lam*X(i))*At + Sig(1)*randn*sqrt(At))*P(i);
        end

        MT(n+1) = M(i+1);
        U(n+1) = log(MT(n));
        
        if MT(n) < 800
            count = count + 1;
        end
        
        hold on
        plot(t,M);
    end
    title('Graph of Mt');
    ylabel('Mt');
    xlabel('t');

    PMt = count/Nmc;
    Exp(j) = sum(U)/Nmc;
    
    if FkExp < Exp(j)
       FkExp = Exp(j);
       MaxTheta = Oo(j);
    end
end

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
figure
plot(Oo,Exp);
title('Graph of Exp');
ylabel('Exp');
xlabel('t');

disp('Probability =');
disp(PMt);

disp('Expectation Value =');
disp(max(Exp));

disp('Theta Value =');
disp(MaxTheta);

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

