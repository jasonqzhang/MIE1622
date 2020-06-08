clear all;
clc
format long;

Nout  = 100000; % number of out-of-sample scenarios
Nin   = 5000;   % number of in-sample scenarios
Ns    = 5;      % number of idiosyncratic scenarios for each systemic

C = 8;          % number of credit states

% Filename to save out-of-sample scenarios
filename_save_out  = 'scen_out';

% Read and parse instrument data
instr_data = dlmread('instrum_data.csv', ',');
instr_id   = instr_data(:,1);           % ID
driver     = instr_data(:,2);           % credit driver
beta       = instr_data(:,3);           % beta (sensitivity to credit driver)
recov_rate = instr_data(:,4);           % expected recovery rate
value      = instr_data(:,5);           % value
prob       = instr_data(:,6:6+C-1);     % credit-state migration probabilities (default to A)
exposure   = instr_data(:,6+C:6+2*C-1); % credit-state migration exposures (default to A)
retn       = instr_data(:,6+2*C);       % market returns

K = size(instr_data, 1); % number of  counterparties

% Read matrix of correlations for credit drivers
rho = dlmread('credit_driver_corr.csv', '\t');
sqrt_rho = (chol(rho))'; % Cholesky decomp of rho (for generating correlated Normal random numbers)


disp('======= Credit Risk Model with Credit-State Migrations =======')
disp('============== Monte Carlo Scenario Generation ===============')
disp(' ')
disp(' ')
disp([' Number of out-of-sample Monte Carlo scenarios = ' int2str(Nout)])
disp([' Number of in-sample Monte Carlo scenarios = ' int2str(Nin)])
disp([' Number of counterparties = ' int2str(K)])
disp(' ')

% Find credit-state for each counterparty
% 8 = AAA, 7 = AA, 6 = A, 5 = BBB, 4 = BB, 3 = B, 2 = CCC, 1 = default
[Ltemp, CS] = max(prob, [], 2);
clear Ltemp

% Account for default recoveries
exposure(:, 1) = (1-recov_rate) .* exposure(:, 1);

% Compute credit-state boundaries
CS_Bdry = norminv( cumsum(prob(:,1:C-1), 2) );

% -------- Insert your code here -------- %

% number of credit drivers is 50
Ncd = length(rho);

if(~exist('scenarios_out.mat','file'))
    
    % -------- Insert your code here -------- 

    % define y with 100,000 systemic scenarios each for 50 correlated credit drivers
    y_MC1 = ones(Nout,Ncd); % size = 100,000 by 50
    % define w as 100,000 systemic scenarios each for 100 counterparties
    w_MC1 = ones(Nout,K); % size = 100,000 by 100
    % define out_of_sample losses as 100,000 systemic scenarios each for
    % the losses of 100 counterparties
    Losses_out = ones(Nout,K); % size = 100,000 by 100
    
    i = 1;
    for s = 1:Nout
        
        y = (randn(1,Ncd) * sqrt_rho)';
        
        % one randomly generated idiosyncratic scenario for each systemic
        z_out = randn(K,1); % size = 100 by 1   
         % for each of the 100 counterparties
             for k = 1:K

                 % calculate creditworthiness index for each scenario and counterparty 
                 w(k) = beta(k) * y(driver(k)) + sqrt(1-beta(k)^2) * z_out(k);
                 % combine w with CS boundaries as a vector in ascending order
                 cs_axis = sort([w(k) CS_Bdry(k,:)]);
                 % returns the position of w amongst the CS boundaries
                 cs_position = find(cs_axis == w(k));
                 % calculate losses for each counterparty based on the position
                 % of w in exposure
                 Losses_out_temp(k) = exposure(k,cs_position);
            
             end
             % store out-of-sample losses, size = 100,000 by 100
             Losses_out(i,:) = Losses_out_temp;
             i = i+1;
    end

    save('scenarios_out', 'Losses_out')
else
    load('scenarios_out', 'Losses_out')
end

% Normal approximation computed from out-of-sample scenarios
mu_l = mean(Losses_out)';
var_l = cov(Losses_out);

% Compute portfolio weights
portf_v = sum(value);     % portfolio value
w0{1} = value / portf_v;  % asset weights (portfolio 1)
w0{2} = ones(K, 1) / K;   % asset weights (portfolio 2)
x0{1} = (portf_v ./ value) .* w0{1};  % asset units (portfolio 1)
x0{2} = (portf_v ./ value) .* w0{2};  % asset units (portfolio 2)

% Quantile levels (99%, 99.9%)
alphas = [0.99 0.999];

% sorting losses for each of the two portfolios
Losses_out_portf{1} = sort(Losses_out * x0{1}); 
Losses_out_portf{2} = sort(Losses_out * x0{2}); 

% Compute VaR and CVaR (non-Normal and Normal) for 100000 scenarios and
% each quantile level
for(portN = 1:2)
    for(q=1:length(alphas))
        alf = alphas(q);
        % -------- Insert your code here -------- %
        VaRout(portN,q)  = Losses_out_portf{portN}(ceil(Nout*alf));
        VaRinN(portN,q)  = mean(Losses_out_portf{portN}) + norminv(alf,0,1)*std(Losses_out_portf{portN});
        CVaRout(portN,q) = (1/(Nout*(1-alf))) * ( (ceil(Nout*alf)-Nout*alf) * VaRout(portN,q) + sum(Losses_out_portf{portN}(ceil(Nout*alf)+1:Nout)));
        CVaRinN(portN,q) = mean(Losses_out_portf{portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(Losses_out_portf{portN});
        % -------- Insert your code here -------- %        
 end
end


% Perform 100 trials
N_trials = 100;    

for(tr=1:N_trials)
    
    % Monte Carlo approximation 1

    % -------- Insert your code here -------- %
    
    i = 1;
    % 1000 systemic scenarios
    for s = 1:ceil(Nin/Ns) 
        
            y_MC1 = (randn(1,Ncd) * sqrt_rho)';
        
            % 5 idiosyncratic scenarios for each systemic
            for si = 1:Ns 
          
            z_MC1 = randn(K,1);  
        
                    % 100 counterparties
                    for k = 1:K
             
                    % calculate creditworthiness index for each scenario and counterparty 
                    w_MC1(k) = beta(k) * y_MC1(driver(k)) + sqrt(1-beta(k)^2) * z_MC1(k);
                    % combine w with CS boundaries as a vector in ascending order
                    cs_axis = sort([w_MC1(k) CS_Bdry(k,:)]);
                    % returns the position of w amongst the CS boundaries
                    cs_position = find(cs_axis == w_MC1(k));
                    % calculate losses for each counterparty based on the position
                    % of w in exposure
                    Losses_inMC1_temp(k) = exposure(k,cs_position);
                    end
              
             % store MC1 losses, size = 5,000 by 100
             Losses_inMC1(i,:) = Losses_inMC1_temp;
             i = i + 1; 
            end
            
    end
        
    
    % Monte Carlo approximation 2
    
    % -------- Insert your code here -------- %
    i = 1;
    % 5000 systemic scenarios (1 idiosyncratic scenario for each systemic)
    for s = 1:Nin 
        
            y_MC2 = (randn(1,Ncd) * sqrt_rho)';
            z_MC2 = randn(K,1);
            
         % for each of the 100 counterparties
         for k = 1:K
             
             % calculate creditworthiness index for each scenario and counterparty 
             w_MC2(k) = beta(k) * y_MC2(driver(k)) + sqrt(1-beta(k)^2) * z_MC2(k);
             % combine w with CS boundaries as a vector in ascending order
             cs_axis = sort([w_MC2(k) CS_Bdry(k,:)]);
             % returns the position of w amongst the CS boundaries
             cs_position = find(cs_axis == w_MC2(k));
             % calculate losses for each counterparty based on the position
             % of w in exposure
             Losses_inMC2_temp(k) = exposure(k,cs_position);
            
         end
             % store MC2 losses, size = 5,000 by 100
             Losses_inMC2(i,:) =  Losses_inMC2_temp;
             i = i+1;
             
    end
      
    
    % Compute VaR and CVaR
    for(portN = 1:2)
        for(q=1:length(alphas))
            alf = alphas(q);
            % -------- Insert your code here -------- %            
            % Compute portfolio loss 
            portf_loss_inMC1{tr,portN} = sort(Losses_inMC1 * x0{portN});
            portf_loss_inMC2{tr,portN} = sort(Losses_inMC2 * x0{portN});
            
            mu_MC1 = mean(Losses_inMC1)';
            var_MC1 = cov(Losses_inMC1);
            mu_MC2 = mean(Losses_inMC2)';
            var_MC2 = cov(Losses_inMC2);
            % Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
            mu_p_MC1 = mu_MC1' * x0{portN};
            sigma_p_MC1 = std(portf_loss_inMC1{portN});
            % Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2
            mu_p_MC2 = mu_MC2' * x0{portN};
            sigma_p_MC2 = std(portf_loss_inMC2{portN});
            % Compute VaR and CVaR for the current trial
            VaRinMC1{portN,q}(tr) = portf_loss_inMC1{portN}(ceil(Nin*alf));
            VaRinMC2{portN,q}(tr) = portf_loss_inMC2{portN}(ceil(Nin*alf));
            VaRinN1{portN,q}(tr) = mean(portf_loss_inMC1{portN}) + norminv(alf,0,1)*std(portf_loss_inMC1{portN});
            VaRinN2{portN,q}(tr) = mean(portf_loss_inMC2{portN}) + norminv(alf,0,1)*std(portf_loss_inMC2{portN});
            CVaRinMC1{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinN1{portN,q}(tr) + sum(portf_loss_inMC1{portN}(ceil(Nin*alf)+1:Nin)));
            CVaRinMC2{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinN2{portN,q}(tr) + sum(portf_loss_inMC2{portN}(ceil(Nin*alf)+1:Nin)));
            CVaRinN1{portN,q}(tr) = mean(portf_loss_inMC1{portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_inMC1{portN});
            CVaRinN2{portN,q}(tr) = mean(portf_loss_inMC2{portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_inMC2{portN});
            % -------- Insert your code here -------- %
        end
    end
end

% Display portfolio VaR and CVaR
for(portN = 1:2)
fprintf('\nPortfolio %d:\n\n', portN)    
 for(q=1:length(alphas))
    alf = alphas(q);
    fprintf('Out-of-sample: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRout(portN,q), 100*alf, CVaRout(portN,q))
    fprintf('In-sample MC1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC1{portN,q}), 100*alf, mean(CVaRinMC1{portN,q}))
    fprintf('In-sample MC2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC2{portN,q}), 100*alf, mean(CVaRinMC2{portN,q}))
    fprintf(' In-sample No: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRinN(portN,q), 100*alf, CVaRinN(portN,q))
    fprintf(' In-sample N1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinN1{portN,q}), 100*alf, mean(CVaRinN1{portN,q}))
    fprintf(' In-sample N2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n\n', 100*alf, mean(VaRinN2{portN,q}), 100*alf, mean(CVaRinN2{portN,q}))
 end
end

% Plot results
%Figure 1
x0 = 10;
y0 = 10;
width = 1000;
height = 400;
% plot portfolio 1 out-of-sample 
figure(1);
set(gcf,'color', 'white','units','points','position',[x0,y0,width,height]);
[frequencyCounts, binLocations] = hist(Losses_out_portf{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
% 99% VaR
line([VaRout(1,1) VaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaR
line([VaRout(1,2) VaRout(1,2)], [0 max(frequencyCounts)/1.5], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaR
line([CVaRout(1,1) CVaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaR
line([CVaRout(1,2) CVaRout(1,2)], [0 max(frequencyCounts)/1.5], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% VaRN
line([VaRinN(1,1) VaRinN(1,1)], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaRN
line([VaRinN(1,2) VaRinN(1,2)], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaRN
line([CVaRinN(1,1) CVaRinN(1,1)], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaRN
line([CVaRinN(1,2) CVaRinN(1,2)], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(Losses_out_portf{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out_portf{1}))/std(Losses_out_portf{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); hold on;

text(0.98*VaRout(1,1), max(frequencyCounts)/1.9, {'VaR','99%'})
text(0.98*VaRout(1,2), max(frequencyCounts)/1.5, {'VaR','99.9%'})
text(0.98*CVaRout(1,1), max(frequencyCounts)/1.9, {'CVaR','99%'})
text(0.98*CVaRout(1,2), max(frequencyCounts)/1.5, {'CVaR','99.9%'})
text(0.93*VaRinN(1,1), max(frequencyCounts)/1.9, {'VaRN','99%'})
text(0.95*VaRinN(1,2), max(frequencyCounts)/1.5, {'VaRN','99.9%'})
text(0.98*CVaRinN(1,1), max(frequencyCounts)/1.9, {'CVaRN','99%'})
text(1*CVaRinN(1,2), max(frequencyCounts)/1.5, {'CVaRN','99.9%'})
hold off;
title('True Distribution of Portfolio 1')
xlabel('Credit-State Migration Exposures')
ylabel('Frequency')


% plot portfolio 2 out-of-sample 
%Figure 2
figure(2);
set(gcf,'color', 'white','units','points','position',[x0,y0,width,height]);
[frequencyCounts, binLocations] = hist(Losses_out_portf{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
% 99% VaR
line([VaRout(2,1) VaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaR
line([VaRout(2,2) VaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaR
line([CVaRout(2,1) CVaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaR
line([CVaRout(2,2) CVaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% VaRN
line([VaRinN(2,1) VaRinN(2,1)], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaRN
line([VaRinN(2,2) VaRinN(2,2)], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaRN
line([CVaRinN(2,1) CVaRinN(2,1)], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaRN
line([CVaRinN(2,2) CVaRinN(2,2)], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(Losses_out_portf{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out_portf{2}))/std(Losses_out_portf{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); hold on;

text(0.98*VaRout(2,1), max(frequencyCounts)/1.9, {'VaR','99%'})
text(0.98*VaRout(2,2), max(frequencyCounts)/1.5, {'VaR','99.9%'})
text(0.98*CVaRout(2,1), max(frequencyCounts)/1.9, {'CVaR','99%'})
text(0.98*CVaRout(2,2), max(frequencyCounts)/1.5, {'CVaR','99.9%'})
text(0.93*VaRinN(2,1), max(frequencyCounts)/1.9, {'VaRN','99%'})
text(0.95*VaRinN(2,2), max(frequencyCounts)/1.5, {'VaRN','99.9%'})
text(0.98*CVaRinN (2,1), max(frequencyCounts)/1.9, {'CVaRN','99%'})
text(1*CVaRinN(2,2), max(frequencyCounts)/1.5, {'CVaRN','99.9%'})
hold off;
title('True Distribution of Portfolio 2')
xlabel('Credit-State Migration Exposures')
ylabel('Frequency')


% plot portfolio 1 in-of-sample 
%Figure 3
VaRinMC1_p1_99 = mean(VaRinMC1{1,1});
VaRinMC1_p1_999 = mean(VaRinMC1{1,2});
CVaRinMC1_p1_99 = mean(CVaRinMC1{1,1});
CVaRinMC1_p1_999 = mean(CVaRinMC1{1,2});
VaRinN1_p1_99 = mean(VaRinN1{1,1});
VaRinN1_p1_999 = mean(VaRinN1{1,2});
CVaRinN1_p1_99 = mean(CVaRinN1{1,1});
CVaRinN1_p1_999 = mean(CVaRinN1{1,2});

% average portfolio loss over 100 trials
total_port_loss_mc1_p1 = zeros(Nin,1);
for i = 1:100

total_port_loss_mc1_p1 = total_port_loss_mc1_p1 + portf_loss_inMC1{i,1};

end

avg_port_loss_mc1_p1 = sort(total_port_loss_mc1_p1/100);   
figure(3);
set(gcf,'color', 'white','units','points','position',[x0,y0,width,height]);
[frequencyCounts, binLocations] = hist(avg_port_loss_mc1_p1, 100);
bar(binLocations, frequencyCounts);
hold on;
% 99% VaR
line([VaRinMC1_p1_99 VaRinMC1_p1_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaR
line([VaRinMC1_p1_999 VaRinMC1_p1_999], [0 max(frequencyCounts)/1.5], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaR
line([CVaRinMC1_p1_99 CVaRinMC1_p1_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaR
line([CVaRinMC1_p1_999 CVaRinMC1_p1_999], [0 max(frequencyCounts)/1.5], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% VaRN
line([VaRinN1_p1_99 VaRinN1_p1_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaRN
line([VaRinN1_p1_999 VaRinN1_p1_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaRN
line([CVaRinN1_p1_99 CVaRinN1_p1_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaRN
line([CVaRinN1_p1_999 CVaRinN1_p1_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(avg_port_loss_mc1_p1)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(avg_port_loss_mc1_p1))/std(avg_port_loss_mc1_p1)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.95*VaRinMC1_p1_99, max(frequencyCounts)/1.9, {'VaR','99%'})
text(0.98*VaRinMC1_p1_999, max(frequencyCounts)/1.5, {'VaR','99.9%'})
text(1*CVaRinMC1_p1_99, max(frequencyCounts)/1.9, {'CVaR','99%'})
text(0.98*CVaRinMC1_p1_999, max(frequencyCounts)/1.5, {'CVaR','99.9%'})
text(0.90*VaRinN1_p1_99, max(frequencyCounts)/1.9, {'VaRN','99%'})
text(0.95*VaRinN1_p1_999, max(frequencyCounts)/1.5, {'VaRN','99.9%'})
text(0.98*CVaRinN1_p1_99, max(frequencyCounts)/1.9, {'CVaRN','99%'})
text(1*CVaRinN1_p1_999, max(frequencyCounts)/1.5, {'CVaRN','99.9%'})
hold off;
title('Distribution of Portfolio 1 (Monte Carlo approximation 1)')
xlabel('Credit-State Migration Exposures')
ylabel('Frequency')

% plot portfolio 2 MC1
%Figure 4
VaRinMC1_p2_99 = mean(VaRinMC1{2,1});
VaRinMC1_p2_999 = mean(VaRinMC1{2,2});
CVaRinMC1_p2_99 = mean(CVaRinMC1{2,1});
CVaRinMC1_p2_999 = mean(CVaRinMC1{2,2});
VaRinN1_p2_99 = mean(VaRinN1{2,1});
VaRinN1_p2_999 = mean(VaRinN1{2,2});
CVaRinN1_p2_99 = mean(CVaRinN1{2,1});
CVaRinN1_p2_999 = mean(CVaRinN1{2,2});

% average portfolio loss over 100 trials
total_port_loss_mc1_p2 = zeros(Nin,1);
for i = 1:100

total_port_loss_mc1_p2 = total_port_loss_mc1_p2 + portf_loss_inMC1{i,1};

end

avg_port_loss_mc1_p2 = sort(total_port_loss_mc1_p2/100);   

figure(4);
set(gcf,'color', 'white','units','points','position',[x0,y0,width,height]);
[frequencyCounts, binLocations] = hist(avg_port_loss_mc1_p2, 100);
bar(binLocations, frequencyCounts);
hold on;
% 99% VaR
line([VaRinMC1_p2_99 VaRinMC1_p2_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaR
line([VaRinMC1_p2_999 VaRinMC1_p2_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaR
line([CVaRinMC1_p2_99 CVaRinMC1_p2_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaR
line([CVaRinMC1_p2_999 CVaRinMC1_p2_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% VaRN
line([VaRinN1_p2_99 VaRinN1_p2_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaRN
line([VaRinN1_p2_999 VaRinN1_p2_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaRN
line([CVaRinN1_p2_99 CVaRinN1_p2_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on
% 99.9% CVaRN
line([CVaRinN1_p2_999 CVaRinN1_p2_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(avg_port_loss_mc1_p2)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(avg_port_loss_mc1_p2))/std(avg_port_loss_mc1_p2)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC1_p2_99, max(frequencyCounts)/1.9, {'VaR','99%'})
text(0.98*VaRinMC1_p2_999, max(frequencyCounts)/1.5, {'VaR','99.9%'})
text(0.98*CVaRinMC1_p2_99, max(frequencyCounts)/1.9, {'CVaR','99%'})
text(0.98*CVaRinMC1_p2_999, max(frequencyCounts)/1.5, {'CVaR','99.9%'})
text(0.95*VaRinN1_p2_99, max(frequencyCounts)/1.9, {'VaRN','99%'})
text(0.95*VaRinN1_p2_999, max(frequencyCounts)/1.5, {'VaRN','99.9%'})
text(1*CVaRinN1_p2_99, max(frequencyCounts)/1.9, {'CVaRN','99%'})
text(1*CVaRinN1_p2_999, max(frequencyCounts)/1.5, {'CVaRN','99.9%'})
hold off;
title('Distribution of Portfolio 2 (Monte Carlo approximation 1)')
xlabel('Credit-State Migration Exposures')
ylabel('Frequency')

% plot portfolio 1 MC2
%Figure 5
VaRinMC2_p1_99 = mean(VaRinMC2{1,1});
VaRinMC2_p1_999 = mean(VaRinMC2{1,2});
CVaRinMC2_p1_99 = mean(CVaRinMC2{1,1});
CVaRinMC2_p1_999 = mean(CVaRinMC2{1,2});
VaRinN2_p1_99 = mean(VaRinN2{1,1});
VaRinN2_p1_999 = mean(VaRinN2{1,2});
CVaRinN2_p1_99 = mean(CVaRinN2{1,1});
CVaRinN2_p1_999 = mean(CVaRinN2{1,2});

% average portfolio loss over 100 trials
total_port_loss_mc2_p1 = zeros(Nin,1);
for i = 1:100

total_port_loss_mc2_p1 = total_port_loss_mc2_p1 + portf_loss_inMC2{i,1};

end

avg_port_loss_mc2_p1 = sort(total_port_loss_mc2_p1/100);  
    
figure(5);
set(gcf,'color', 'white','units','points','position',[x0,y0,width,height]);
[frequencyCounts, binLocations] = hist(avg_port_loss_mc2_p1, 100);
bar(binLocations, frequencyCounts);
hold on;
% 99% VaR
line([VaRinMC2_p1_99 VaRinMC2_p1_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaR
line([VaRinMC2_p1_999 VaRinMC2_p1_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaR
line([CVaRinMC2_p1_99 CVaRinMC2_p1_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaR
line([CVaRinMC2_p1_999 CVaRinMC2_p1_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% VaRN
line([VaRinN2_p1_99 VaRinN2_p1_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaRN
line([VaRinN2_p1_999 VaRinN2_p1_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaRN
line([CVaRinN2_p1_99 CVaRinN2_p1_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on
% 99.9% CVaRN
line([CVaRinN2_p1_999 CVaRinN2_p1_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;


normf = ( 1/(std(avg_port_loss_mc2_p1)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(avg_port_loss_mc2_p1))/std(avg_port_loss_mc2_p1)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC2_p1_99, max(frequencyCounts)/1.9, {'VaR','99%'})
text(0.98*VaRinMC2_p1_999, max(frequencyCounts)/1.5, {'VaR','99.9%'})
text(0.98*CVaRinMC2_p1_99, max(frequencyCounts)/1.9, {'CVaR','99%'})
text(0.98*CVaRinMC2_p1_999, max(frequencyCounts)/1.5, {'CVaR','99.9%'})
text(0.95*VaRinN2_p1_99, max(frequencyCounts)/1.9, {'VaRN','99%'})
text(0.95*VaRinN2_p1_999, max(frequencyCounts)/1.5, {'VaRN','99.9%'})
text(1*CVaRinN2_p1_99, max(frequencyCounts)/1.9, {'CVaRN','99%'})
text(1*CVaRinN2_p1_999, max(frequencyCounts)/1.5, {'CVaRN','99.9%'})
hold off;
title('Distribution of Portfolio 1 (Monte Carlo approximation 2)')
xlabel('Credit-State Migration Exposures')
ylabel('Frequency')

% plot portfolio 2 MC2
%Figure 6
VaRinMC2_p2_99 = mean(VaRinMC2{2,1});
VaRinMC2_p2_999 = mean(VaRinMC2{2,2});
CVaRinMC2_p2_99 = mean(CVaRinMC2{2,1});
CVaRinMC2_p2_999 = mean(CVaRinMC2{2,2});
VaRinN2_p2_99 = mean(VaRinN2{2,1});
VaRinN2_p2_999 = mean(VaRinN2{2,2});
CVaRinN2_p2_99 = mean(CVaRinN2{2,1});
CVaRinN2_p2_999 = mean(CVaRinN2{2,2});

% average portfolio loss over 100 trials
total_port_loss_mc2_p2 = zeros(Nin,1);
for i = 1:100

total_port_loss_mc2_p2 = total_port_loss_mc2_p2 + portf_loss_inMC2{i,1};

end

avg_port_loss_mc2_p2 = sort(total_port_loss_mc2_p2/100);  

figure(6);
% set(gcf, 'color', 'white');
set(gcf,'color', 'white','units','points','position',[x0,y0,width,height])
[frequencyCounts, binLocations] = hist(avg_port_loss_mc2_p2, 100);
bar(binLocations, frequencyCounts);
hold on;
% 99% VaR
line([VaRinMC2_p2_99 VaRinMC2_p2_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaR
line([VaRinMC2_p2_999 VaRinMC2_p2_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaR
line([CVaRinMC2_p2_99 CVaRinMC2_p2_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% CVaR
line([CVaRinMC2_p2_999 CVaRinMC2_p2_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% VaRN
line([VaRinN2_p2_99 VaRinN2_p2_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99.9% VaRN
line([VaRinN2_p2_999 VaRinN2_p2_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;
% 99% CVaRN
line([CVaRinN2_p2_99 CVaRinN2_p2_99], [0 max(frequencyCounts)/2], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '--');
hold on
% 99.9% CVaRN
line([CVaRinN2_p2_999 CVaRinN2_p2_999], [0 max(frequencyCounts)/1.5], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(avg_port_loss_mc2_p2)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(avg_port_loss_mc2_p2))/std(avg_port_loss_mc2_p2)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC2_p2_99, max(frequencyCounts)/1.9, {'VaR','99%'})
text(0.98*VaRinMC2_p2_999, max(frequencyCounts)/1.5, {'VaR','99.9%'})
text(0.98*CVaRinMC2_p2_99, max(frequencyCounts)/1.9, {'CVaR','99%'})
text(0.98*CVaRinMC2_p2_999, max(frequencyCounts)/1.5, {'CVaR','99.9%'})
text(0.95*VaRinN2_p2_99, max(frequencyCounts)/1.9, {'VaRN','99%'})
text(0.95*VaRinN2_p2_999, max(frequencyCounts)/1.5, {'VaRN','99.9%'})
text(1*CVaRinN2_p2_99, max(frequencyCounts)/1.9, {'CVaRN','99%'})
text(1*CVaRinN2_p2_999, max(frequencyCounts)/1.5, {'CVaRN','99.9%'})
hold off;
title('Distribution of Portfolio 2 (Monte Carlo approximation 2)')
xlabel('Credit-State Migration Exposures')
ylabel('Frequency')