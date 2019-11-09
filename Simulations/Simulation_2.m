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
At = T/N;

%Initialization of coordinates
W(1) = 0; %Brownian Motion
S(1) = 10; %for underlying asset
t = (0:N)*At;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%         Underlying Asset Monte Carlo Simulation     %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

for k = 1:Nmc
    for i = 1:N
        g = randn();
        W(i+1) = W(i) + g * sqrt(At);
        S(i+1) = S(i) * ((1+(r*At)) + sigma*g*sqrt(At));
    end
    FinalValueMB(k) = W(end);
end
expectation = sum(FinalValueMB)/Nmc;
disp('Expectation ='); disp(expectation);