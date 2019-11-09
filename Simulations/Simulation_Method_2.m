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
sigma = 0.5;
N = 100;
Nmc = 1000;
At = 0.005;

%Initialization of coordinates
for k=1:Nmc
    W(1,k) = 0; %Brownian Motion
    S(1,k) = 10; %for underlying asset
end
t = (0:N)*At;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%         Underlying Asset Monte Carlo Simulation     %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

for k = 1:Nmc
    for i = 1:N
        g = randn;
        S(i+1,k) = S(i,k) * exp((r-((sigma^2)/2))*At + sqrt(At)*sigma*g);
    end
end

X = S(N+1,:);

figure
plot(t,S);
title('t => S');
ylabel('St');
xlabel('t');

stat = Stats(X);
disp('Expectation of Asset =');
disp(stat(1));

disp('Variance of Asset =');
disp(stat(2));

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    FUNCTIONS                        %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%Function parameters
function [EV] = Stats(x)
    %Variables
    n = length(x);
    
    %Core of function
    Exp = sum(x)/n;
    Var = sum((x - Exp).^2)/n;
    
    %The results
    EV = [Exp;Var];
    
end
