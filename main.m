% Undergraduate Thesis -- Xinyi Hu
% 50ETF Implied Volatility  
% Department of Finance, School of Economics, Peking University

close all; clear; clc;

%% Plot Settings
set(groot, 'DefaultLineLineWidth', 2.5, ...
           'DefaultTextInterpreter', 'LaTeX', ...
           'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
           'DefaultAxesFontName', 'LaTeX', ...
           'DefaultLegendInterpreter', 'LaTeX', ...
           'DefaultAxesLineWidth', 1.5, ...
           'DefaultAxesFontSize', 16, ...
           'DefaultAxesBox', 'on', ...
           'DefaultAxesColor', [1, 1, 1], ...
           'DefaultFigureColor', [1, 1, 1], ...
           'DefaultFigureInvertHardcopy', 'off', ...
           'DefaultFigurePaperUnits', 'inches', ...
           'DefaultFigureUnits', 'inches', ...
           'DefaultFigurePaperPosition', [0, 0, 5, 3.09], ...
           'DefaultFigurePaperSize', [10, 6], ...
           'DefaultFigurePosition', [0.1, 0.1, 4.9, 2.99]);

%% Table Combanition

DailyPrice = readtable('./data/Option/DailyPrice.xlsx');
Info = readtable('./data/Option/Info.xlsx');
Shibor = readtable('./data/InterestRate/Shibor.xlsx');
ETF = readtable('./data/50ETF/510050.xlsx');
TradingCalendar = readtable('./data/50ETF/TradingDays.xlsx');

% Convert shibor actual/360 compounding to continuous compunding
% Save as risk free rate table Rf
Rf = table();
Rf.Date = Shibor.Date;
Rf.Rf_1D = log(1 + 0.01*Shibor.S_1D*1/360) * 360/1;
Rf.Rf_1W = log(1 + 0.01*Shibor.S_1W*7/360) * 360/7;
Rf.Rf_2W = log(1 + 0.01*Shibor.S_2W*14/360) * 360/14;
Rf.Rf_1M = log(1 + 0.01*Shibor.S_1M*30/360) * 360/30;
Rf.Rf_3M = log(1 + 0.01*Shibor.S_3M*90/360) * 360/90;
Rf.Rf_6M = log(1 + 0.01*Shibor.S_6M*180/360) * 360/180;
Rf.Rf_9M = log(1 + 0.01*Shibor.S_9M*270/360) * 360/270;
Rf.Rf_1Y = log(1 + 0.01*Shibor.S_1Y*360/360) * 360/360;

% Include Info
Full_1 = innerjoin(DailyPrice,Info,'Keys',{'Code','ExecCode','Name'});
% Include ETF close price
Full_2 = innerjoin(Full_1,ETF,'Keys',{'Date'},'RightVariables',{'ETFClose'});
% Include Rf
Full = innerjoin(Full_2,Rf,'Keys',{'Date'});
% Sort by date and option name
Full = sortrows(Full,{'Date','Name'},{'ascend','ascend'});

clear Full_1 Full_2;

%% Variable Construction

% Calculate days to maturity (and in year), only trading days are counted
FullTradingDays = unique(Full.Date);
Full.DaysToMaturity = days_to_maturity(Full.Date, Full.MatureDate, TradingCalendar.Date);
Full.YearToMaturity = Full.DaysToMaturity / 250;

% Calculate risk free rate w.r.t DaysToMaturity (linear interpolation)
Full.Rf = rf_shibor(Full.YearToMaturity, Full.Rf_1D, Full.Rf_1W, Full.Rf_2W,...
    Full.Rf_1M, Full.Rf_3M, Full.Rf_6M, Full.Rf_9M, Full.Rf_1Y);

% Calculate dividend adjusted ETFClose
Full.D = zeros(size(Full.Date,1),1);
Full.D(cellfun(@(s)contains(s, 'A'), Full.ExecCode)) = 0.05;
Full.ETFClosePrime = Full.ETFClose - Full.D .* exp(-Full.Rf .* Full.YearToMaturity);
Full.q = zeros(size(Full.Date,1),1);

% Calculate the moneyness (Strike price / Forward price)
Full.ForwardPrice = exp(Full.Rf .* Full.YearToMaturity) .* Full.ETFClosePrime;
Full.Moneyness = Full.Strike ./ Full.ForwardPrice;

% Calculate minimum and maximum theoretical price
Full.CallMin = Full.ETFClosePrime - ...
    Full.Strike .* exp(-Full.Rf .* Full.YearToMaturity);
Full.PutMin = Full.Strike .* exp(-Full.Rf .* Full.YearToMaturity) - ...
    Full.ETFClosePrime;
Full.CallMax = Full.ETFClosePrime;
Full.PutMax = Full.Strike .* exp(-Full.Rf .* Full.YearToMaturity);


%% Observation Filtering

% Visualize the volume & amount by moneyness & YearToMaturity
%Call = Full(strcmp(Full.Type,'EuropeanCall'),:);
%Put = Full(strcmp(Full.Type,'EuropeanPut'),:);
subplot(2,2,1);
scatter(Full.Moneyness, Full.Volume);
title('50ETF Option Volume by Moneyness');
subplot(2,2,2);
scatter(Full.Moneyness, Full.Amount);
title('50ETF Option Amount by Moneyness');
subplot(2,2,3);
scatter(Full.YearToMaturity, Full.Volume);
title('50ETF Option Volume by YearToMaturity');
subplot(2,2,4);
scatter(Full.YearToMaturity, Full.Amount);
title('50ETF Option Amount by YearToMaturity');

% Description statistics Volume/Amount vs Moneyness
thres = [0.0,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,10];
for i=1:length(thres)-1
    fprintf('Moneyness:%4.1f<M<=%4.1f  Obs:%d  Mean Volume:%8.2f  Mean Amount:%12.2f\n',...
    thres(i), thres(i+1),...
    length(Full.Volume(thres(i) < Full.Moneyness & Full.Moneyness <= thres(i+1))),...
    mean(Full.Volume(thres(i) < Full.Moneyness & Full.Moneyness <= thres(i+1)),'omitnan'),...
    mean(Full.Amount(thres(i) < Full.Moneyness & Full.Moneyness <= thres(i+1)),'omitnan'));
end
fprintf('All Obs Mean Volume: %8.2f\n',mean(Full.Volume,'omitnan'));
fprintf('All Obs Mean Amount: %8.2f\n',mean(Full.Amount,'omitnan'));
clear thres;

% Drop options with moneyness < 0.7 or > 1.5
% Drop ITM options
% Drop options with YearToMaturity <= 5/250 (a week)
% Drop options out of price bound
% Drop options with volume < 100
% Drop options close 
f1 = Full.Moneyness >= 0.7 & Full.Moneyness <= 1.5;
f2 = (Full.Moneyness >= 1 & strcmp(Full.Type,'EuropeanCall')) |...
    (Full.Moneyness < 1 & strcmp(Full.Type,'EuropeanPut'));
f3 = Full.YearToMaturity > 5/250;
f4 = (strcmp(Full.Type,'EuropeanCall') & Full.Close >= Full.CallMin) |...
    (strcmp(Full.Type,'EuropeanPut') & Full.Close >= Full.PutMin);
f5 = (strcmp(Full.Type,'EuropeanCall') & Full.Close <= Full.CallMax) |...
    (strcmp(Full.Type,'EuropeanPut') & Full.Close <= Full.PutMax);
f6 = Full.Volume >= 50;
f7 = (Full.Close < Full.HighLimit) & (Full.Close > Full.LowLimit);
Filtered = Full(f1 & f2 & f3 & f4 & f5 & f6 & f7, :);
FilteredTradingDays = unique(Filtered.Date);

clear f1 f2 f3 f4 f5 f6 f7;

% Plot num of Obs 
Obs = zeros(size(FullTradingDays,1),2);
for i=1:size(FullTradingDays,1)
    Obs(i,1) = sum(Full.Date == FullTradingDays(i));
    Obs(i,2) = sum(Filtered.Date == FullTradingDays(i));
end
plot(FullTradingDays, Obs(:,1), 'color','k','LineWidth',1);
hold on
plot(FullTradingDays, Obs(:,2),'--', 'color','k','LineWidth',0.6);
legend('Before Filtration','After Filtration')
%title('Number of Observations over Time')
hold off

%% Market Implied Volatility

Filtered.MarketIV = bsmiv(Filtered.Type, Filtered.Close, Filtered.ETFClosePrime,...
    Filtered.Strike, Filtered.YearToMaturity, Filtered.Rf, Filtered.q);

% Calculate volume weighted BSMIV of call/put to yield the volatility
% surface
[G,jDate,jMoneyness,jYearToMaturity] = findgroups(...
    Filtered.Date,Filtered.Moneyness,Filtered.YearToMaturity);
MarketIVSurface = table();
MarketIVSurface.Date = jDate; 
MarketIVSurface.Moneyness = jMoneyness;
MarketIVSurface.YearToMaturity = jYearToMaturity;
MarketIVSurface.IV = splitapply(@weightediv,Filtered(:,{'MarketIV','Volume'}),G);
MarketIVSurface.LogIV = log(MarketIVSurface.IV);

clear G jDate jMoneyness jYearToMaturity;

% Some IVS(discrete) plots 
for i=60:60:960
    tmp = MarketIVSurface(MarketIVSurface.Date == FilteredTradingDays(i), :);
    subplot(4,4,i/60);
    %[sf, gof]  = fit([x, y],z,'poly23');
    scatter3(tmp.Moneyness,tmp.YearToMaturity,tmp.IV);
    title(string(FilteredTradingDays(i)));
    xlabel('Moneyness','FontSize',12); ylabel('YTM','FontSize',12); zlabel('MarketIV','FontSize',12);
    %set(gca,'XLim',[0.7 1.5],'YLim',[0 1],'ZLim',[0 0.5]);
end

% Summary stats of MarketIV
tmp = Filtered(Filtered.Moneyness < 0.95 & Filtered.YearToMaturity < 0.1,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness < 0.95 & Filtered.YearToMaturity >= 0.1 & Filtered.YearToMaturity <= 0.3,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness < 0.95 & Filtered.YearToMaturity > 0.3,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness >= 0.95 & Filtered.Moneyness <= 1.05 & Filtered.YearToMaturity < 0.1,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness >= 0.95 & Filtered.Moneyness <= 1.05 & Filtered.YearToMaturity >= 0.1 & Filtered.YearToMaturity <= 0.3,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness >= 0.95 & Filtered.Moneyness <= 1.05 & Filtered.YearToMaturity > 0.3,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness > 1.05 & Filtered.YearToMaturity < 0.1,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness > 1.05 & Filtered.YearToMaturity >= 0.1 & Filtered.YearToMaturity <= 0.3,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));
tmp = Filtered(Filtered.Moneyness > 1.05 & Filtered.YearToMaturity > 0.3,:).MarketIV;
fprintf('obs:%5d  mean:%4.4f  std:%4.4f\n',length(tmp),mean(tmp),std(tmp));


%% Implied Volatility Modelling (Heston)

% Calibrate Heston params
HestonWindow = 3;

HestonParams = table('Size',[size(FilteredTradingDays,1)+1-HestonWindow,6],...
    'VariableTypes',["datetime","double","double","double","double","double"],...
    'VariableNames',["Date","v0","theta","rho","kappa","sigma"]);
HestonParams.Date = FilteredTradingDays(LinearWindow:end);

startparams = [0.2^2, 2, -0.5, 0.25^2, 0.3];

for i=1:size(FilteredTradingDays,1)+1-HestonWindow
    tmp = Filtered(Filtered.Date >= FilteredTradingDays(i) &...
        Filtered.Date <= FilteredTradingDays(i+HestonWindow-1),:);
    res = lsqnonlin(@(x) hestonprice(...
        tmp.Type,tmp.ETFClosePrime,tmp.Strike,tmp.YearToMaturity,tmp.Rf,tmp.q,x(1),x(2),x(3),x(4),x(5)) - tmp.Close,...
        startparams, ...% v0,theta,rho,kappa,sigma
        [eps eps -1+eps eps eps  ], ... % lower bound for parameter vector
        [Inf Inf 1-eps Inf  Inf  ], ...  % upper bound for parameter vector
        optimoptions('lsqnonlin'));
    HestonParams(i,2:end) = num2cell(res);
end

clear startparmas tmp res;

% Calculate Heston IV
HestonFiltered = innerjoin(Filtered(:,{'Date','Type','ETFClosePrime','Strike','YearToMaturity','Rf','q','Close','MarketIV'}), ...
    HestonParams, 'Keys',{'Date'});
HestonFiltered.HestonPrice = cell2mat(rowfun(@hestonprice, ...
    HestonFiltered(:,{'Type','ETFClosePrime','Strike','YearToMaturity','Rf','q','v0','theta','rho','kappa','sigma'}),...
    'OutputFormat','cell'));
HestonFiltered.HestonIV = bsmiv(HestonFiltered.Type, HestonFiltered.HestonPrice,...
    HestonFiltered.ETFClosePrime, HestonFiltered.Strike, ...
    HestonFiltered.YearToMaturity, HestonFiltered.Rf, HestonFiltered.q);
HestonFiltered.residual = HestonFiltered.MarketIV - HestonFiltered.HestonIV;
[G]=findgroups(HestonFiltered.Date);
HestonParams.rmse = splitapply(@(x)sqrt(mean(x.^2)), HestonFiltered.residual, G);
HestonParams.rmse_avg = HestonParams.rmse ./ splitapply(@mean, HestonFiltered.MarketIV, G);

plot(FilteredTradingDays(HestonWindow:end),movmean(HestonParams.rmse_avg, 1),'k','LineWidth',1);
ylabel('RMSE/AVG');
ylim([0,0.8]);


%% Implied Volatility Modelling (Linear)

LinearWindow = 3;

% Moneyness order 1 YearToMaturity order 1
LinearParams_11 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,7],...
    'VariableTypes',["datetime","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","adjrsquare","rmse","rmse_avg"]);
LinearParams_11.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly11');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_11(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 1 YearToMaturity order 2
LinearParams_12 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,9],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P11","P02","adjrsquare","rmse","rmse_avg"]);
LinearParams_12.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly12');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_12(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p11,sf.p02,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 2 YearToMaturity order 1
LinearParams_21 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,9],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P20","P11","adjrsquare","rmse","rmse_avg"]);
LinearParams_21.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly21');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_21(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 2 YearToMaturity order 2
LinearParams_22 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,10],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P20","P11","P02","adjrsquare","rmse","rmse_avg"]);
LinearParams_22.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly22');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_22(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,sf.p02,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 1 YearToMaturity order 3
LinearParams_13 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,11],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P11","P02","P12","P03","adjrsquare","rmse","rmse_avg"]);
LinearParams_13.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly13');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_13(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p11,sf.p02,sf.p12,sf.p03,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 3 YearToMaturity order 1
LinearParams_31 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,11],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P20","P11","P30","P21","adjrsquare","rmse","rmse_avg"]);
LinearParams_31.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly31');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_31(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,sf.p30,sf.p21,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 2 YearToMaturity order 3
LinearParams_23 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,13],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P20","P11","P02","P21","P12","P03","adjrsquare","rmse","rmse_avg"]);
LinearParams_23.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly23');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_23(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,sf.p02,sf.p21,sf.p12,sf.p03,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 3 YearToMaturity order 2
LinearParams_32 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,13],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P20","P11","P02","P30","P21","P12","adjrsquare","rmse","rmse_avg"]);
LinearParams_32.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly32');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_32(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,sf.p02,sf.p30,sf.p21,sf.p12,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

% Moneyness order 3 YearToMaturity order 3
LinearParams_33 = table('Size',[size(FilteredTradingDays,1)+1-LinearWindow,14],...
    'VariableTypes',["datetime","double","double","double","double","double","double","double","double","double","double","double","double","double"],...
    'VariableNames',["Date","P00","P10","P01","P20","P11","P02","P30","P21","P12","P03","adjrsquare","rmse","rmse_avg"]);
LinearParams_33.Date = FilteredTradingDays(LinearWindow:end);

for i=1:size(FilteredTradingDays,1)+1-LinearWindow
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly33');
    residual = tmp.IV - exp(sf([tmp.Moneyness, tmp.YearToMaturity]));
    LinearParams_33(i,2:end) = num2cell([sf.p00,sf.p10,sf.p01,sf.p20,sf.p11,sf.p02,sf.p30,sf.p21,sf.p12,sf.p03,gof.adjrsquare,...
        sqrt(mean(residual .^ 2, 'omitnan')),...
        sqrt(mean(residual .^ 2, 'omitnan')) / mean(tmp.IV, 'omitnan')]);
    %plot(sf,[tmp.Moneyness, tmp.YearToMaturity],tmp.MarketIV);
end

clear indx tmp sf gof residual;

plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_11.adjrsquare,60),'k','LineWidth',1);
hold on 
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_12.adjrsquare,60),'k-.','LineWidth',1);
hold on 
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_21.adjrsquare,60),'k:','LineWidth',2);
hold on
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_22.adjrsquare,60),'k','LineWidth',2);
hold on
legend({'Poly11','Poly12','Ploy21','Poly22'});
ylabel('Adjusted R-Square Moving Average');
ylim([0,1]);
hold off

plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_13.adjrsquare,60),'k','LineWidth',1);
hold on
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_31.adjrsquare,60),'k--','LineWidth',1);
hold on
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_23.adjrsquare,60),'k:','LineWidth',2);
hold on
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_32.adjrsquare,60),'k','LineWidth',2);
hold on
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_33.adjrsquare,60),'k-.','LineWidth',1);
legend({'Poly13','Poly31','Poly23','Poly32','Poly33'});
ylabel('Adjusted R-Square Moving Average');
ylim([0,1]);
hold off

plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_11.adjrsquare,60),'k:','LineWidth',2);
hold on 
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_21.adjrsquare,60),'k-.','LineWidth',1);
hold on 
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_31.adjrsquare,60),'k','LineWidth',2);
hold on
legend({'Poly11','Poly21','Poly31'});
ylabel('Adjusted R-Square Moving Average');
ylim([0,1]);
hold off

mean(LinearParams_21.adjrsquare);
std(LinearParams_21.adjrsquare);
plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_21.adjrsquare,1),'k','LineWidth',2);

%% Error Comparison

plot(FilteredTradingDays(LinearWindow:end),movmean(LinearParams_21.rmse_avg, 30), 'k:','LineWidth',2);
hold on
plot(FilteredTradingDays(HestonWindow:end),movmean(HestonParams.rmse_avg, 30),'k','LineWidth',1);
ylabel('RMSE/AVG Moving Average');
ylim([0,0.5]);
legend({"Linear -- Poly21","Heston Model"});
hold off

% Some IVS(continuous) plots 

m = 0.7:0.01:1.5;
t = 0.02:0.01:0.67;
[M,T] = meshgrid(m,t);

% Moneyness order 2 YearToMaturity order 1
for i=60+1-LinearWindow:60:960+1-LinearWindow
    subplot(4,4,(i+LinearWindow-1)/60);
    indx = (MarketIVSurface.Date >= FilteredTradingDays(i)) & ...
        (MarketIVSurface.Date <= FilteredTradingDays(i+LinearWindow-1));
    tmp = MarketIVSurface(indx, :);
    [sf, gof]  = fit([tmp.Moneyness, tmp.YearToMaturity],tmp.LogIV,'poly21');
    surf(M,T,exp(sf.p00+sf.p10.*M+sf.p01.*T+sf.p20.*M.^2+sf.p11.*M.*T));
    title(string(FilteredTradingDays(i+LinearWindow-1)));
    xlabel('Moneyness','FontSize',12); ylabel('YTM','FontSize',12); zlabel('MarketIV','FontSize',12);
end

clear indx tmp sf gof;

% Heston
for i=60+1-HestonWindow:60:960+1-HestonWindow
    subplot(4,4,(i+LinearWindow-1)/60);
    tmp = Filtered(Filtered.Date == FilteredTradingDays(i+HestonWindow-1),:);
    res = table2array(tmp(1,{'ETFClosePrime','Rf','q'}));
    ETFClosePrime=res(1); rf=res(2); q=res(3); 
    res=table2array(HestonParams(i,2:end));
    v0=res(1);theta=res(2);rho=res(3);kappa=res(4);sigma=res(5);
    K = ETFClosePrime * exp(rf.*T) .* M; 
    iv = nan(size(K,1),size(K,2));
    for j=1:size(K,1)
        for k=1:size(K,2)
            if M(j,k)>=1
                type = "EuropeanCall";
            else
                type = "EuropeanPut";
            end
            price = hestonprice(type,ETFClosePrime,K(j,k),T(j,k),rf,q,v0,theta,rho,kappa,sigma);
            if isnan(price)
                disp("Heston price not calculated")
            elseif price < 0.0001
                flag = 1;
            else
                iv(j,k) = bsmiv(type,price,ETFClosePrime,K(j,k),T(j,k),rf,q);
            end
        end
    end
    surf(M,T,iv);
    title(string(FilteredTradingDays(i+LinearWindow-1)));
    xlabel('Moneyness','FontSize',12); ylabel('YTM','FontSize',12); zlabel('MarketIV','FontSize',12);
end

clear tmp res ETFClosePrime v0 sigma theta rho kappa K iv rf q type price;
clear m t M T;

%% Param Plot

% plot poly21 params
subplot(5,1,1);
plot(LinearParams_21.Date, LinearParams_21.P00, 'k','LineWidth',1);
legend('P00');
ylim([-80,80]);
subplot(5,1,2);
plot(LinearParams_21.Date, LinearParams_21.P10, 'k','LineWidth',1);
legend('P10');
ylim([-100,100]);
subplot(5,1,3);
plot(LinearParams_21.Date, LinearParams_21.P01, 'k','LineWidth',1);
legend('P01');
ylim([-20,20]);
subplot(5,1,4);
plot(LinearParams_21.Date, LinearParams_21.P20, 'k','LineWidth',1);
legend('P20');
ylim([-80,80]);
subplot(5,1,5);
plot(LinearParams_21.Date, LinearParams_21.P11, 'k','LineWidth',1);
legend('P11');
ylim([-20,20]);

% plot poly21 params
subplot(5,1,1);
plot(HestonParams.Date, HestonParams.v0, 'k','LineWidth',1);
legend('v0');
ylim([0,0.5]);
subplot(5,1,2);
plot(HestonParams.Date, HestonParams.theta, 'k','LineWidth',1);
legend('theta');
ylim([0,5]);
subplot(5,1,3);
plot(HestonParams.Date, HestonParams.rho, 'k','LineWidth',1);
legend('rho');
ylim([-2,2]);
subplot(5,1,4);
plot(HestonParams.Date, HestonParams.kappa, 'k','LineWidth',1);
legend('kappa');
ylim([0,20]);
subplot(5,1,5);
plot(HestonParams.Date, HestonParams.sigma, 'k','LineWidth',1);
legend('sigma');
ylim([0,5]);


%% Parmas Predicability

% Sample and out-sample

s = 800;

% VAR Estimate
% Poly21 in sample
numseries = 5;
LinearParamArray = table2array(LinearParams_21(1:s,{'P00','P10','P01','P20','P11'}));
for p=1:10
    Mdl = varm(numseries, p);
    EstMdl = estimate(Mdl, LinearParamArray);
    res = summarize(EstMdl);
    disp(res.BIC);
end

p=2;
LinearMdl = varm(numseries, p);
[LinearEstMdl,LinearEstSE,LinearlogL,LinearInsampleResidual] = estimate(LinearMdl, LinearParamArray);
%summarize(EstMdl);
LinearInsampleError = abs(LinearInsampleResidual);
LinearInsampleRelativeError = LinearInsampleError ./ abs(LinearParamArray(3:800,:));
NaiveInsampleError = abs(LinearParamArray(2:799,:) - LinearParamArray(3:800,:));
NaiveInsampleRelativeError = NaiveInsampleError ./ abs(LinearParamArray(3:800,:));
mean(LinearInsampleError)
mean(NaiveInsampleError)

% Poly21 out of sample
c = LinearEstMdl.Constant;
A1 = cell2mat(LinearEstMdl.AR(1));
A2 = cell2mat(LinearEstMdl.AR(2)); 
LinearParamArrayO = table2array(LinearParams_21(s-1:end,{'P00','P10','P01','P20','P11'}));
LinearOutsampleMV = c' + LinearParamArrayO(1:181,:) * A1' + LinearParamArrayO(2:182,:) * A2';
LinearOutsampleError = abs(LinearOutsampleMV - LinearParamArrayO(3:end,:));
LinearOutsampleRelativeError = LinearOutsampleError ./ abs(LinearParamArrayO(3:end,:));

NaiveOutsampleError = abs(LinearParamArrayO(2:182,:) - LinearParamArrayO(3:end,:));
NaiveOutsampleRelativeError = NaiveOutsampleError ./ abs(LinearParamArrayO(3:end,:));

mean(LinearOutsampleError)
mean(NaiveOutsampleError)


% Heston in sample
numseries = 5;
HestonParamArray = table2array(HestonParams(1:s,{'v0','theta','rho','kappa','sigma'}));
for p=1:10
    Mdl = varm(numseries, p);
    EstMdl = estimate(Mdl, HestonParamArray);
    res = summarize(EstMdl);
    disp(res.BIC);
end

p=7;
HestonMdl = varm(numseries, p);
[HestonEstMdl,HestonEstSE,HestonlogL,HestonInsampleResidual] = estimate(HestonMdl, HestonParamArray);
%summarize(EstMdl);
HestonInsampleError = abs(HestonInsampleResidual);
HestonInsampleRelativeError = HestonInsampleError ./ abs(HestonParamArray(8:800,:));
NaiveInsampleError = abs(HestonParamArray(7:799,:) - HestonParamArray(8:800,:));
NaiveInsampleRelativeError = NaiveInsampleError ./ abs(HestonParamArray(8:800,:));
mean(HestonInsampleError)
mean(NaiveInsampleError)

% Poly21 out of sample
c = HestonEstMdl.Constant;
A1 = cell2mat(HestonEstMdl.AR(1));
A2 = cell2mat(HestonEstMdl.AR(2)); 
A3 = cell2mat(HestonEstMdl.AR(3));
A4 = cell2mat(HestonEstMdl.AR(4)); 
A5 = cell2mat(HestonEstMdl.AR(5));
A6 = cell2mat(HestonEstMdl.AR(6)); 
A7 = cell2mat(HestonEstMdl.AR(7)); 

HestonParamArrayO = table2array(HestonParams(s-6:end,{'v0','theta','rho','kappa','sigma'}));
HestonOutsampleMV = c' + HestonParamArrayO(1:181,:) * A1' + HestonParamArrayO(2:182,:) * A2' + HestonParamArrayO(3:183,:) * A3' + HestonParamArrayO(4:184,:) * A4' + HestonParamArrayO(5:185,:) * A5' + HestonParamArrayO(6:186,:) * A6' + HestonParamArrayO(7:187,:) * A7';
HestonOutsampleError = abs(HestonOutsampleMV - HestonParamArrayO(8:end,:));
HestonOutsampleRelativeError = HestonOutsampleError ./ abs(HestonParamArrayO(8:end,:));

NaiveOutsampleError = abs(HestonParamArrayO(7:187,:) - HestonParamArrayO(8:end,:));
NaiveOutsampleRelativeError = NaiveOutsampleError ./ abs(HestonParamArrayO(8:end,:));

mean(HestonOutsampleError)
mean(NaiveOutsampleError)
