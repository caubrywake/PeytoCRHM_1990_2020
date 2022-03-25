%% A2 - streamflow benchmark
%% compare measured streamflow to internannual avarge 
 % our model has to be better than the simplest possible model which is the average of all the years
 
close all
clear all

addpath(genpath('D:\PeytoCRHM_1990_2020\chrm'))
addpath('D:\PeytoCRHM_1990_2020\function')
savedir = 'D:\PeytoCRHM_1990_2020\data_process\modeleval\'
figdir = 'D:\PeytoCRHM_1990_2020\fig\modeleval\'

load('D:\SaltDilution\F_buildingratingcurve\output\Q_SR50_daily_20132020.mat')
bf_meas = Q_SR50_daily;
t_meas = Q_SR50_daily_time+days(1);
clear Q_SR50_daily Q_SR50_daily_time


%% remove leap days for measured
% measured flow
vectime= datevec(t_meas);
leapdays = find(vectime(:, 2) == 2 & vectime(:, 3) == 29);
bf_meas(leapdays)= [];
t_meas(leapdays) = [];
clear vectime leapdays

%% Reshape streamflow for annual values
% add nan to have 2020 complete
t_meas_pad = [t_meas; (t_meas(end)+days(1):days(1):datetime('15-Oct-2020 00:00:00'))'];
bf_meas_pad = [bf_meas; nan(length(t_meas(end)+days(1):days(1):datetime('15-Oct-2020')), 1)];

% average for 15 may-15 oct
vectime =  datevec(t_meas_pad);
idx_may15 = find(vectime(:, 2)== 05 & vectime(:, 3) == 15);
idx_oct15 = find(vectime(:, 2)== 10 & vectime(:, 3) == 15);

t_meas_pad(idx_may15) % just to check
t_meas_pad(idx_oct15) % just to check

for i = 1:length(idx_may15) % not the last year
bfmeas_doy(:, i) = bf_meas_pad(idx_may15(i):idx_oct15(i))
end 

bfmeas_doy(bfmeas_doy==0)=nan; % change 0 to nan
bfmeas_doy_mean = nanmean(bfmeas_doy, 2);

figure;
p1  = plot(t_meas(idx_may15(1):idx_oct15(1)), bfmeas_doy, 'k');
hold on
p2 = plot(t_meas(idx_may15(1):idx_oct15(1)), bfmeas_doy_mean, 'r', 'linewidth', 1.2);
xlabel ('DOY');
ylabel('Measured Streamflow (m3/s)')
legend ([p1(1) p2(1)], 'Individual years (2013-2020)', 'Average')

%% Create synthetic benchmarck
% where May 15-Oct 15 every year is the interannual doy average
bf_bench = bf_meas_pad;

for i = 1:length(idx_may15) % not the last year
bf_bench(idx_may15(i):idx_oct15(i)) = bfmeas_doy_mean;
end 

plot(t_meas_pad, bf_bench)
%% Set benchmark as the modelled results
t_bench = t_meas_pad;
bf_bench = bf_bench;
%% Join together to get 
% select dates with flow
flowdates  = [datetime('04-Jun-2013') datetime('1-Oct-2013');... 
             datetime('7-June-2014') datetime('14-Oct-2014');...
             datetime('22-May-2015') datetime('16-Oct-2015');...
             datetime('14-May-2016') datetime('10-Oct-2016');...
             datetime('27-May-2017') datetime('13-Oct-2017');...
             datetime('18-May-2018') datetime('01-Oct-2018');...
             datetime('27-May-2019') datetime('7-Oct-2019');...
             datetime('4-Jun-2020') datetime('17-Sep-2020')];
 flowdates_column = reshape(flowdates, 16, 1);
 % extract the correposnding dates for the measured and benchelled flow
 for i = 1:length(flowdates_column)
idx_flowdates_meas(i) = find(t_meas == flowdates(i));
idx_flowdates_bench(i) = find(t_bench == flowdates(i));
 end 

 idx_flowdates_meas = reshape(idx_flowdates_meas, 8, 2);
 idx_flowdates_bench = reshape(idx_flowdates_meas, 8, 2);
 % compile the average of all the DOY
 
  %% Calculate benchel evaluation statistic for each year
 clear bench meas tbench tmeas RMSEyr NSEyr MByr MAEyr R2yr
 for i = 1:length(flowdates)
 bench = bf_bench(idx_flowdates_bench(i, 1): idx_flowdates_bench(i, 2));
 meas = bf_meas( idx_flowdates_meas(i, 1): idx_flowdates_meas(i, 2));
 
 tbench = datenum(t_bench( idx_flowdates_bench(i, 1): idx_flowdates_bench(i, 2)));
 tmeas = datenum(t_meas( idx_flowdates_meas(i, 1): idx_flowdates_meas(i, 2)));
 % rmse
 RMSEyrbench(i, 1) =  sqrt(mean((bench - meas).^2));
 % r2
 r = corrcoef(bench,meas);  R2yrbench(i, 1) = r(2)^2;
 % Mean absloute error
 MAEyrbench(i, 1) = mean(abs(bench - meas));
 %Nash sutcliffe
[NSEyrbench(i, 1), metricid] = nashsutcliffe([tmeas meas], [tbench bench]);
% mean Bias
MByrbench(i, 1) = (sum(bench-meas)/sum(meas));
 end 
 
 % Make a table
yr = (2013:2020)';
Q_Stats_perYear_benchmark = table (yr, NSEyrbench, RMSEyrbench, R2yrbench, MAEyrbench, MByrbench);
headers = {'year', 'NSE', 'RMSE', 'R2', 'MAE', 'MB'}; 
Q_Stats_perYear.Properties.VariableNames = headers;
Q_Stats_perYear_benchmark 
writetable(Q_Stats_perYear_benchmark, strcat(savedir,'Q_Stats_perYear_benchmark', datestr(now, 'yyyymmddHHMM')))

%% For the entire period
% keep only times with flow
bench = bf_bench;
meas = bf_meas;
tbench = t_bench;
tmeas = t_meas;

figure; 
plot(tmeas, meas); hold on; plot(tbench, bench)

for i = 1:length(idx_flowdates_meas)-1 % change all the values outside of measured period as nan
   bench(idx_flowdates_bench(i, 2):idx_flowdates_bench(i+1, 1)) = nan;
   meas(idx_flowdates_meas(i, 2):idx_flowdates_meas(i+1, 1)) = nan;
end 

% cut to just 2013-2020
bench = bench(idx_flowdates_bench(1,1):idx_flowdates_bench(8,2));
tbench = tbench(idx_flowdates_bench(1,1):idx_flowdates_bench(8,2));
meas = meas(idx_flowdates_meas(1,1):idx_flowdates_meas(8,2));
tmeas = tmeas(idx_flowdates_meas(1,1):idx_flowdates_meas(8,2));
% remove all the nan
a = find(isnan(bench));
tbench(a)=[];
bench(a)=[];
a = find(isnan(meas));
tmeas(isnan(meas)) = [];
meas(isnan(meas))=[];

figure; 
plot(tmeas, meas); hold on; plot(tbench, bench)
% Calculate statistics
RMSEbench =  sqrt(mean((bench - meas).^2)); % rmse
r = corrcoef(bench,meas);  R2bench = r(2)^2; % r2
MAEbench = mean(abs(bench - meas)); % Mean absloute error
[NSEbench, metricid] = nashsutcliffe([datenum(tmeas) meas], [datenum(tbench) bench]); %Nash sutcliffe
MBbench = (sum(bench-meas)/sum(meas));% mean Bias

% add ti previous table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Figures;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% All year scatter pot
figure; 
sc1 = scatter(meas, bench, '.k'); hold on
ls = lsline;
rf = refline(1,0)
ls(1).Color = 'r';
rf(1).Color = 'k';
xlabel ('Daily Mean Streamflow, Measured (m^3/s)');
ylabel('Daily Mean Streamflow, benchelled (m^3/s)');
legend('obs vs bench', 'least square line','1:1 line', 'location', 'Southeast') ;
grid on
text (.2,11, {strcat('R2=',num2str(R2bench));strcat('RMSE=',num2str(RMSEbench));strcat('NSE=',num2str(NSEbench))});

figname ='Basinflow_scatterplot_benchmark';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))
%%
fig = figure('units','inches','outerposition',[0 0 6 7]);
subplot(4,2,1)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
ylabel ({'streamflow, m^3/s'})
grid on
ylim([0 10])
xlim ([datetime('01-May-2013') datetime('15-Oct-2013')])
text (datetime('2-May-2013'),9, {strcat('NSE=',num2str(round(NSEyrbench(1),2)))});
 legend ( 'Meas', 'bench')

subplot(4,2,2)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2014') datetime('15-Oct-2014')])
text (datetime('2-May-2014'),9, {strcat('NSE=',num2str(round(NSEyrbench(2),2)))});

subplot(4,2,3)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
grid on
xlim ([datetime('01-May-2015') datetime('15-Oct-2015')])
ylabel ({'streamflow, m^3/s'})
ylim([0 10])
text (datetime('2-May-2015'),9, {strcat('NSE=',num2str(round(NSEyrbench(3),2)))})

subplot(4,2,4)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
ylabel ({'streamflow, m^3/s'})
grid on
ylim([0 10])
xlim ([datetime('01-May-2016') datetime('15-Oct-2016')])
text (datetime('2-May-2016'),9, {strcat('NSE=',num2str(round(NSEyrbench(4),2)))});

subplot(4,2,5)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2017') datetime('15-Oct-2017')])
text (datetime('2-May-2017'),9, {strcat('NSE=',num2str(round(NSEyrbench(5),2)))});

subplot(4,2,6)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
grid on
xlim ([datetime('01-May-2018') datetime('15-Oct-2018')])
ylabel ({'streamflow, m^3/s'})
ylim([0 10])
text (datetime('2-May-2018'),9, {strcat('NSE=',num2str(round(NSEyrbench(6),2)))})

subplot(4,2,7)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
ylabel ({'streamflow, m^3/s'})
grid on
% legend ( 'Measured', 'benchelled','Orientation', 'Horizontal', 'location', 'NorthOutside')
ylim([0 10])
xlim ([datetime('01-May-2019') datetime('15-Oct-2019')])
text (datetime('2-May-2019'),9, {strcat('NSE=',num2str(round(NSEyrbench(7),2)))});

subplot(4,2,8)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_bench , bf_bench, 'r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2020') datetime('15-Oct-2020')])
text (datetime('2-May-2020'),9, {strcat('NSE=',num2str(round(NSEyrbench(8),2)))});

figname ='Basinflow_perYr_benchmark';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))

