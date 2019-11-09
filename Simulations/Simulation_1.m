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
At = T/N;

%Initialization of coordinates
W(1) = 0; %for Brownian Motion
S(1) = 10; %for underlying asset
t = (0:N)*At;

%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                    MAIN PROGRAM                     %
%         Underlying Asset Monte Carlo Simulation     %
%/////////////////////////-\\\\\\\\\\\\\\\\\\\\\\\\\\\\

for i = 1:N
    g = randn;
    W(i+1) = W(i) + g * sqrt(At);
    S(i+1) = S(i) * ((1+(r*At)) + sigma*g*sqrt(At));
end
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

figure
plot(t,W,t,S);
title('t => W & t => S');
ylabel('St');
xlabel('t');