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
        W(i+1,k) = W(i,k) + g * sqrt(At);
        S(i+1,k) = S(i,k) * ((1+(r*At)) + sigma*g*sqrt(At));
    end
end
figure
plot(W);
title('W');
ylabel('St');
xlabel('t');

figure
plot(t,W);
title('t => W');
ylabel('St');
xlabel('t');

figure
plot(t,S);
title('t => S');
ylabel('St');
xlabel('t');