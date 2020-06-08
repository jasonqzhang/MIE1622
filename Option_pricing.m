clc;
clear all;
format long

% Pricing a European option using Black-Scholes formula and Monte Carlo simulations
% Pricing a Barrier option using Monte Carlo simulations

S0 = 100;     % spot price of the underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.1;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
Sb = 110;     % barrier

% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations
numPaths = 50000;
numSteps = 12;

% Implement your Black-Scholes pricing formula
[call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);

% Implement your one-step Monte Carlo pricing procedure for European option
% numSteps = 1;
[callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for European option
[callMC_European_Price_multi_step, putMC_European_Price_multi_step, MCpaths] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths);

% Implement your one-step Monte Carlo pricing procedure for Barrier option
% numSteps = 1;
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);

disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(putBS_European_Price)])
disp(['One-step MC price of an European call option is ',num2str(callMC_European_Price_1_step)])
disp(['One-step MC price of an European put option is ',num2str(putMC_European_Price_1_step)])
disp(['Multi-step MC price of an European call option is ',num2str(callMC_European_Price_multi_step)])
disp(['Multi-step MC price of an European put option is ',num2str(putMC_European_Price_multi_step)])
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of  an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])


% Plot results
figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
for i=1:numPaths 
    plot(1:numSteps+1,MCpaths(:,i));
    hold on;
end
hold off;
axis([0 numSteps+1 30 inf])
title('Underlying Stock Price Simulations');

%% Black-Scholes Equation 
function [call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma)

t = 0;
d1 = (1/sigma*sqrt(T-t)) * (log(S0/K) + (r+sigma^2/2)*(T-t));
d2 = d1 - sigma*sqrt(T-t);
call_BS_European_Price = normcdf(d1)*S0 - normcdf(d2)*K*exp(-r*(T-t));
putBS_European_Price = normcdf(-d2)*K*exp(-r*(T-t)) - normcdf(-d1)*S0;

end

%% Multi-Step MC
function [callMC_European_Price_multi_step, putMC_European_Price_multi_step, paths] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths)

    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    paths(1,:) = S0;

        % Generates paths with the corresponding steps
        for iPath = 1:numPaths
            for iStep = 1:numSteps
                paths(iStep+1, iPath) = paths(iStep, iPath) * exp((mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1));
            end
        end 

        call = zeros(numPaths,1);
        put = zeros(numPaths,1);
        
        % Calculates Call and Put prices for each path
        for iPath = 1:numPaths
         call(iPath,1) = max(paths(numSteps+1,iPath) - K, 0) * exp(-r*T);
         put(iPath,1) = max(K - paths(numSteps+1,iPath), 0) * exp(-r*T);
        end

     % Price of the option is the average of all paths
     callMC_European_Price_multi_step = mean(call);
     putMC_European_Price_multi_step = mean(put);
    
end

%% Multi-Step MC Barrier Knock-in
% The option becomes a standard option if the barrier price was crossed
% somtime before expiration T

function [callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)

    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    paths(1,:) = S0;
    
    % Generates paths with the corresponding steps
     for iPath = 1:numPaths
            for iStep = 1:numSteps
                paths(iStep+1, iPath) = paths(iStep, iPath) * exp((mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1));
            end
     end 
        
      call = zeros(numPaths,1);
      put = zeros(numPaths,1);
      
     % Calculates Call and Put prices for each path based on the barrier
     % Option will be valid if price at any point before expiration
     % crosses the barrier price Sb
     for iPath = 1:numPaths
         if any(paths(:,iPath) > Sb)
             
             call(iPath,1) = max(paths(numSteps+1,iPath) - K, 0) * exp(-r*T);
             put(iPath,1) = max(K - paths(numSteps+1,iPath), 0) * exp(-r*T);
              
         else
             call(iPath,1) = 0;
             put(iPath,1) = 0;
         end
         
     end
     
    % Price of the option is the average of all paths 
    callMC_Barrier_Knockin_Price_multi_step = mean(call);
    putMC_Barrier_Knockin_Price_multi_step = mean(put);
    
end

