clc;
clear all;
format long

global Q cur_year init_positions 
% Input files
% 2015-2016
input_file_prices  = 'Daily_closing_prices.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);

% Find the number of trading days in Nov-Dec 2014 and
% compute expected return and covariance matrix for period 1
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
cur_returns0 = data_prices(day_ind_start0+1:day_ind_end0,:) ./ data_prices(day_ind_start0:day_ind_end0-1,:) - 1;
mu = mean(cur_returns0)';
Q = cov(cur_returns0);

% Remove datapoints for year 2014
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Initial positions in the portfolio
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';

% Initial value of the portfolio
init_value = data_prices(1,:) * init_positions;
fprintf('\nInitial portfolio value = $ %10.2f\n\n', init_value);

% Initial portfolio weights
w_init = (data_prices(1,:) .* init_positions')' / init_value;

% Number of periods, assets, trading days
N_periods = 6*length(unique(dates_array(:,1))); % 6 periods per year
N = length(tickers);
N_days = length(dates);

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;
% Annual risk-free rate for years 2008-2009 is 4.5%
r_rf2008_2009 = 0.045;

% Number of strategies
strategy_functions = {'strat_buy_and_hold' 'strat_equally_weighted' 'strat_min_variance' 'strat_max_Sharpe' 'strat_equal_risk_contr' 'strat_lever_equal_risk_contr' 'strat_robust_optim'};
strategy_names     = {'Buy and Hold' 'Equally Weighted Portfolio' 'Minimum Variance Portfolio' 'Maximum Sharpe Ratio Portfolio' 'Equal Risk Contributions Portfolio' 'Leveraged Equal Risk Contributions Portfolio' 'Robust Optimization Portfolio'};
%N_strat = 5; % comment this in your code
N_strat = length(strategy_functions); % uncomment this in your code
fh_array = cellfun(@str2func, strategy_functions, 'UniformOutput', false);

for (period = 1:N_periods)
   % Compute current year and month, first and last day of the period
   if(dates_array(1,1)==15)
       cur_year  = 15 + floor(period/7);
   else
       cur_year  = 2015 + floor(period/7);
   end
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   fprintf('\nPeriod %d: start date %s, end date %s\n', period, char(dates(day_ind_start)), char(dates(day_ind_end)));

   % Prices for the current day
   cur_prices = data_prices(day_ind_start,:);

   % Execute portfolio selection strategies
   for(strategy = 1:N_strat)

      % Get current portfolio positions
      if(period==1)
         curr_positions = init_positions;
         curr_cash = 0;
         portf_value{strategy} = zeros(N_days,1);
      else
         curr_positions = x{strategy,period-1};
         curr_cash = cash{strategy,period-1};
      end

      % Compute strategy
      [x{strategy,period} cash{strategy,period} weight{strategy,period}] = fh_array{strategy}(curr_positions, curr_cash, mu, Q, cur_prices);

      if strategy == 6
          
          % Compute portfolio value for strategy 6 excluding borrowed assets
          portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period} - init_value;
      else
          % Compute portfolio value for other strategies
          portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period};
      end 
      fprintf('   Strategy "%s", value begin = $ %10.2f, value end = $ %10.2f\n', char(strategy_names{strategy}), portf_value{strategy}(day_ind_start), portf_value{strategy}(day_ind_end));

   end
      
   % Compute expected returns and covariances for the next period
   cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
   mu = mean(cur_returns)';
   Q = cov(cur_returns);
   
end

% Plotting results

% Graph 1: daily portfolio value of each strat. for 504 trading days
figure(1)
% Note that interest payments are already included in the leverage strategy
plot(1:N_days,portf_value{1,1},1:N_days,portf_value{1,2},1:N_days,portf_value{1,3},1:N_days,portf_value{1,4},1:N_days,portf_value{1,5},1:N_days,portf_value{1,6},1:N_days,portf_value{1,7});
title('Performance of Strategies for 2015 and 2016');
xlabel('Trading Days'); ylabel('Portfolio Value'); legend('Buy and Hold','Equally Weighted Portfolio','Minimum Variance Portfolio','Maximum Sharpe Ratio Portfolio','Equal Risk Contributions Portfolio','Leveraged Equal Risk Contributions Portfolio','Robust Optimization Portfolio');

% Graph 2: dynamic changes in portfolio allocations for strat_min_variance
figure(2)
% The weight function output is converted into a 20x13 matrix for strat_3 
weight_3 = [w_init weight{3,1} weight{3,2} weight{3,3} weight{3,4} weight{3,5} weight{3,6} weight{3,7} weight{3,8} weight{3,9} weight{3,10} weight{3,11} weight{3,12}];
 
% Each row(stock) of weight_3 is plotted against the number of periods
        for i = 1:N
         plot(0:N_periods, weight_3(i,:));
         hold on
        end 
        
ylim([0 1])
title('Min Variance Strategy Dynamic Change'); xlabel('Period'); ylabel('Weight'); legend('1','2','3','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20');

% Graph 3: dynamic changes in portfolio allocations for strat_max_Sharpe
figure(3)
% The weight function output is converted into a 20x13 matrix for strat_4
weight_4 = [w_init weight{4,1} weight{4,2} weight{4,3} weight{4,4} weight{4,5} weight{4,6} weight{4,7} weight{4,8} weight{4,9} weight{4,10} weight{4,11} weight{4,12}];
% Each row(stock) of weight_4 is plotted against the number of periods
        for i = 1:N
         plot(0:N_periods, weight_4(i,:));
         hold on
        end 
        
ylim([0 1])
title('Max Sharpe Ratio Strategy Dynamic Change'); xlabel('Period'); ylabel('Weight'); legend('1','2','3','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20');

% Graph 4: dynamic changes in portfolio allocations for strat_robust_optim
figure(4)
% The weight function output is converted into a 20x13 matrix for strat_7 
weight_7 = [w_init weight{7,1} weight{7,2} weight{7,3} weight{7,4} weight{7,5} weight{7,6} weight{7,7} weight{7,8} weight{7,9} weight{7,10} weight{7,11} weight{7,12}];
% Each row(stock) of weight_7 is plotted against the number of periods
        for i = 1:N
         plot(0:N_periods, weight_7(i,:));
         hold on
        end 
ylim([0 1])
title('Robust Mean-Variance Strategy Dynamic Change'); xlabel('Period'); ylabel('Weight'); legend('1','2','3','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20');