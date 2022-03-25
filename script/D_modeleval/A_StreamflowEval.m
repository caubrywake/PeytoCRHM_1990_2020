%% 1987-2020 model evaluatioln

%% set-up
close all
clear all

addpath(genpath('D:\4_PeytoCRHM_1990_2020\chrm'))
addpath('D:\4_PeytoCRHM_1990_2020\function')
savedir = 'D:\4_PeytoCRHM_1990_2020\data_process\modeleval\'
figdir = 'D:\4_PeytoCRHM_1990_2020\fig\modeleval\'

load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\basinflow.mat')
 load('D:\1_SaltDilution\E_buildingratingcurve\output\Q_SR50_daily_2013_2020.mat')
bf_meas = Q_SR50_daily;
t_meas = Q_SR50_daily_time+days(1);
clear Q_SR50_daily Q_SR50_daily_time


%% Change modelled flow to daily
T = timetable(timeCRHM, basinflow/3600, basingw/3600);
TT = retime (T, 'daily', 'mean');
bf_mod = table2array(TT);
t_mod = TT.timeCRHM;

%% Quick visual comparison
figure
plot(t_mod, bf_mod(:, 1));hold on
plot(t_mod, bf_mod(:, 2));
plot(t_meas, bf_meas); 
legend ('modelled sw', 'mod gw',  'measured sw')
ylabel ('Daily Streamflow (m3/s)')

%% remove leap days for measured and modelled
% modeleld flow
vectime= datevec(t_mod);
leapdays = find(vectime(:, 2) == 2 & vectime(:, 3) == 29);
t_mod(leapdays,:)= [];
bf_mod(leapdays, :)= [];
% measured flow
vectime= datevec(t_meas);
leapdays = find(vectime(:, 2) == 2 & vectime(:, 3) == 29);
bf_meas(leapdays)= [];
t_meas(leapdays) = [];
clear vectime leapdays

%% Reshape streamflow for annual values

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
 % extract the correposnding dates for the measured and modelled flow
 for i = 1:length(flowdates_column)
 idx_flowdates_meas(i) = find(t_meas == flowdates(i));
 idx_flowdates_mod(i) = find(t_mod == flowdates(i));
 end 

 idx_flowdates_meas = reshape(idx_flowdates_meas, 8, 2);
 idx_flowdates_mod = reshape(idx_flowdates_mod, 8, 2);
 %% Calculate Model evaluation statistic for each year
 gw_mod = bf_mod(:, 2);
 bf_mod = bf_mod(:, 1);
 

 clear mod meas tmod tmeas RMSEyr NSEyr MByr MAEyr R2yr
 for i = 1:length(flowdates)
 mod = bf_mod(idx_flowdates_mod(i, 1): idx_flowdates_mod(i, 2));
 meas = bf_meas( idx_flowdates_meas(i, 1): idx_flowdates_meas(i, 2));
 
 tmod = datenum(t_mod( idx_flowdates_mod(i, 1): idx_flowdates_mod(i, 2)));
 tmeas = datenum(t_meas( idx_flowdates_meas(i, 1): idx_flowdates_meas(i, 2)));
 % rmse
 RMSEyr(i, 1) =  sqrt(mean((mod - meas).^2));
 KGyr(i, 1) = kinggupta(meas, mod);
 % r2
 r = corrcoef(mod,meas);  R2yr(i, 1) = r(2)^2;
 % Mean absloute error
 MAEyr(i, 1) = mean(abs(mod - meas));
 %Nash sutcliffe
[NSEyr(i, 1), metricid] = nashsutcliffe([tmeas meas], [tmod mod]);
% mean Bias
MByr(i, 1) = (sum(mod-meas)/sum(meas));
 end 
 
 % Make a table
yr = (2013:2020)';
Q_Stats_perYear = table (yr, NSEyr, RMSEyr, R2yr, MAEyr, MByr, KGyr);
headers = {'year', 'NSE', 'RMSE', 'R2', 'MAE', 'MB', 'KG'}; 
Q_Stats_perYear .Properties.VariableNames = headers;
Q_Stats_perYear 
writetable(Q_Stats_perYear, strcat(savedir,'Q_Stats_perYear_', datestr(now, 'yyyymmddHHMM')))

%% For the entire period
% keep only times with flow
mod = bf_mod;
meas = bf_meas;
tmod = t_mod;
tmeas = t_meas;

for i = 1:length(idx_flowdates_meas)-1 % change all the values outside of measured period as nan
   mod(idx_flowdates_mod(i, 2):idx_flowdates_mod(i+1, 1)) = nan;
   meas(idx_flowdates_meas(i, 2):idx_flowdates_meas(i+1, 1)) = nan;
end 

% cut to just 2013-2020
mod = mod(idx_flowdates_mod(1,1):idx_flowdates_mod(8,2));
tmod = tmod(idx_flowdates_mod(1,1):idx_flowdates_mod(8,2));
meas = meas(idx_flowdates_meas(1,1):idx_flowdates_meas(8,2));
tmeas = tmeas(idx_flowdates_meas(1,1):idx_flowdates_meas(8,2));
% remove all the nan
a = find(isnan(mod));
tmod(a)=[];
mod(a)=[];
a = find(isnan(meas));
tmeas(isnan(meas)) = [];
meas(isnan(meas))=[];

% Calculate statistics
RMSE =  sqrt(mean((mod - meas).^2)); % rmse
r = corrcoef(mod,meas);  R2 = r(2)^2; % r2
MAE = mean(abs(mod - meas)); % Mean absloute error
[NSE, metricid] = nashsutcliffe([datenum(tmeas) meas], [datenum(tmod) mod]); %Nash sutcliffe
MB = (sum(mod-meas)/sum(meas));% mean Bias
KG = kinggupta(meas, mod);
% add ti previous table
yr = {'2013';'2014';'2015';'2016'; '2017';'2018';'2019';'2020' ; 'all'};
Q_Stats = table (yr, [NSEyr; NSE], [RMSEyr ;RMSE], [R2yr ;R2], [MAEyr ;MAE], [MByr; MB], [KGyr;KG]);
headers = {'year', 'NSE', 'RMSE', 'R2', 'MAE', 'MB', 'KG'}; 
Q_Stats.Properties.VariableNames = headers;
writetable(Q_Stats, strcat(savedir,'Q_Stats_', datestr(now, 'yyyymmddHHMM')))
Q_Stats 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Figures;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% All year scatter pot
fig = figure('units','inches','outerposition',[0 0 4 4]);
fig.Renderer='Painters';
sc1 = scatter(meas, mod, '.k'); hold on
ls = lsline;
rf = refline(1,0)
ls(1).Color = 'r';
rf(1).Color = 'k';
R = corrcoef(meas, mod, 'rows', 'pairwise'); R2 = r(2) * r(2);
xlabel ('Daily Mean Streamflow, Measured (m^3/s)');
ylabel('Daily Mean Streamflow, Modelled (m^3/s)');
legend([ls(1) rf(1)], 'ls line','1:1 line', 'location', 'Southeast') ;
grid on
xlim([0 12]); ylim([0 12]);
text (.2,10, {strcat('r^2=',num2str(round(R2,2)));strcat('RMSE=',num2str(round(RMSE,2)));strcat('NSE=',num2str(round(NSE,2)))});

figname ='Basinflow_scatterplot';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))
%%
fig = figure('units','inches','outerposition',[0 0 8 9]);
fig.Renderer='Painters';
sp1= subplot(4,2,1)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod, gw_mod, ':r')
ylabel ({'Q (m^3/s)'})
grid on
ylim([0 10])
xlim ([datetime('01-May-2013') datetime('1-Mar-2014')])
% xtickformat('dd-MMM')
text (datetime('2-May-2013'),9, strcat('(a)',{' '}, ...
    'NSE=',num2str(round(NSEyr(1),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(1), 2)), ',' , {' '},'MB = ', num2str(round(MByr(1),2))))
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

lg =  legend ( 'Meas. Streamflow', 'Mod. Streamflow', 'Mod. Groundwater','orientation', 'vertical','location', 'northeast')
pos = get(lg,'Position');
lg.Position = [0.2635    0.8218    0.2004    0.0785];

sp1 = subplot(4,2,2)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2014') datetime('1-Mar-2015')])
text (datetime('2-May-2014'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(2),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(2), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(2),2))))
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

sp1=subplot(4,2,3)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
grid on
xlim ([datetime('01-May-2015') datetime('1-Mar-2016')])
ylabel ({'Q (m^3/s)'})
ylim([0 10])
text (datetime('2-May-2015'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(3),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(3), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(3),2))));
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

sp1 = subplot(4,2,4)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2016') datetime('1-Mar-2017')])
text (datetime('2-May-2016'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(4),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(4), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(4),2))));
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

sp1 = subplot(4,2,5)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2017') datetime('1-Mar-2018')])
text (datetime('2-May-2017'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(5),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(5), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(5),2))));
ylabel ({'Q (m^3/s)'})
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

sp1 = subplot(4,2,6)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
grid on
xlim ([datetime('01-May-2018') datetime('1-Mar-2019')])
ylim([0 10])
text (datetime('2-May-2018'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(6),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(6), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(6),2))));
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

sp1 = subplot(4,2,7)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
ylabel ({'Q (m^3/s)'})
grid on
ylim([0 10])
xlim ([datetime('01-May-2019') datetime('1-Mar-2020')])
text (datetime('2-May-2019'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(7),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(7), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(7),2))));
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

sp1 = subplot(4,2,8)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2020') datetime('1-Mar-2021')])
text (datetime('2-May-2020'),9,strcat( ...
    'NSE=',num2str(round(NSEyr(8),2)), ',' , {' '},...
    'KG =', num2str(round(KGyr(8), 2)), ',' , {' '},...
    'MB = ', num2str(round(MByr(8),2))));
possb = get(sp1,'Position');
set(sp1,'Position',[possb(1)-0.05 possb(2) possb(3)+0.05 possb(4)+0.02]);

figname ='Basinflow_perYr';
saveas (gcf, strcat(figdir, figname,  datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))

%% as one continuous
fig = figure('units','inches','outerposition',[0 0 10 5]);
bf_meas_bckup=bf_meas;
bf_meas_bckup(bf_meas_bckup<=0)=nan;
plot(t_meas, bf_meas_bckup, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
ylabel ({'Q (m^3/s)'})
grid on
ylim([0 10.3])
 legend ( 'Measured', 'Modelled- streamflow', 'Modelled - groundwater flow', 'location', 'northoutside','orientation', 'horizontal')
xlim ([datetime('01-May-2013') datetime('15-Oct-2020')])
text (datetime('2-May-2013'),6.5, {strcat('NSE=',num2str(round(NSE,2)))});

%% Calculate benchmark
% bf_meas = bf_meas_bckup;
% a = find(isnan(bf_meas))
% bf_meas(a)==0;
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

%% Set benchmark as the modelled results
t_bench = t_meas_pad;
bf_bench = bf_bench;

 for i = 1:length(flowdates_column)
idx_flowdates_meas(i) = find(t_meas == flowdates(i));
idx_flowdates_bench(i) = find(t_bench == flowdates(i));
 end 

 idx_flowdates_bench = reshape(idx_flowdates_meas, 8, 2);
 
  %% Calculate benchel evaluation statistic for each year
 for i = 1:length(flowdates)
 bench = bf_bench(idx_flowdates_bench(i, 1): idx_flowdates_bench(i, 2));
 meas = bf_meas( idx_flowdates_meas(i, 1): idx_flowdates_meas(i, 2));
%  meas(isnan(meas))=0;
 
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
% Calculate statistics
RMSEbench =  sqrt(mean((bench - meas).^2)); % rmse
r = corrcoef(bench,meas);  R2bench = r(2)^2; % r2
MAEbench = mean(abs(bench - meas)); % Mean absloute error
[NSEbench, metricid] = nashsutcliffe([datenum(tmeas) meas], [datenum(tbench) bench]); %Nash sutcliffe
MBbench = (sum(bench-meas)/sum(meas));% mean Bias
%% Grah with benchmark
%%
fig = figure('units','inches','outerposition',[0 0 6 7]);
fig.Renderer='Painters';
subplot(4,2,1)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
ylabel ({'Q (m^3/s)'})
grid on
ylim([0 10])
xlim ([datetime('01-May-2013') datetime('15-Oct-2013')])
text (datetime('2-May-2013'),9, {strcat('NSE=',num2str(round(NSEyr(1),2)))});
text (datetime('2-May-2013'),7, {strcat('NSEb=',num2str(round(NSEyrbench(1),2)))});
 legend ( 'Meas', 'Mod', 'Benchmark')

subplot(4,2,2)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2014') datetime('15-Dec-2014')])
text (datetime('2-May-2014'),9, {strcat('NSE=',num2str(round(NSEyr(2),2)))});
text (datetime('2-May-2014'),7, {strcat('NSEb=',num2str(round(NSEyrbench(2),2)))});

subplot(4,2,3)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
grid on
xlim ([datetime('01-May-2015') datetime('15-Oct-2015')])
ylabel ({'Q (m^3/s)'})
ylim([0 10])
text (datetime('2-May-2015'),9, {strcat('NSE=',num2str(round(NSEyr(3),2)))})
text (datetime('2-May-2015'),7, {strcat('NSEb=',num2str(round(NSEyrbench(3),2)))});

subplot(4,2,4)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
ylabel ({'streamflow, m^3/s'})
grid on
ylim([0 10])
xlim ([datetime('01-May-2016') datetime('15-Oct-2016')])
text (datetime('2-May-2016'),9, {strcat('NSE=',num2str(round(NSEyr(4),2)))});
text (datetime('2-May-2016'),7, {strcat('NSEb=',num2str(round(NSEyrbench(4),2)))});

subplot(4,2,5)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2017') datetime('15-Oct-2017')])
text (datetime('2-May-2017'),9, {strcat('NSE=',num2str(round(NSEyr(5),2)))});
text (datetime('2-May-2017'),7, {strcat('NSEb=',num2str(round(NSEyrbench(5),2)))});
ylabel ({'Q (m^3/s)'})

subplot(4,2,6)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
grid on
xlim ([datetime('01-May-2018') datetime('15-Oct-2018')])
ylim([0 10])
text (datetime('2-May-2018'),9, {strcat('NSE=',num2str(round(NSEyr(6),2)))})
text (datetime('2-May-2018'),7, {strcat('NSEb=',num2str(round(NSEyrbench(6),2)))});

subplot(4,2,7)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
ylabel ({'Q (m^3/s)'})
grid on
% legend ( 'Measured', 'Modelled','Orientation', 'Horizontal', 'location', 'NorthOutside')
ylim([0 10])
xlim ([datetime('01-May-2019') datetime('15-Oct-2019')])
text (datetime('2-May-2019'),9, {strcat('NSE=',num2str(round(NSEyr(7),2)))});
text (datetime('2-May-2019'),7, {strcat('NSEb=',num2str(round(NSEyrbench(7),2)))});

subplot(4,2,8)
plot(t_meas, bf_meas, 'k'); hold on
plot(t_mod , bf_mod, 'r'); hold on
plot(t_bench , bf_bench, 'Color', [.5 .5 .5]); hold on
grid on
ylim([0 10])
xlim ([datetime('01-May-2020') datetime('15-Oct-2020')])
text (datetime('2-May-2020'),9, {strcat('NSE=',num2str(round(NSEyr(8),2)))});
text (datetime('2-May-2020'),7, {strcat('NSEb=',num2str(round(NSEyrbench(8),2)))});

figname ='Basinflow_perYr_withbench';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))

%% as one continuous
fig = figure('units','inches','outerposition',[0 0 10 3]);
fig.Renderer='Painters';
plot(t_meas, smooth(bf_meas, 10), 'k'); hold on
plot(t_mod , smooth(bf_mod, 10), 'r'); hold on
plot(t_mod , gw_mod, ':r'); hold on
ylabel ({'Q (m^3/s)'})
grid on
ylim([0 7])
 legend ( 'Meas', 'Mod- SW', 'Mod - GW', 'location', 'northoutside','orientation', 'horizontal')
xlim ([datetime('01-May-2013') datetime('15-Oct-2020')])
text (datetime('2-May-2013'),6.5, {strcat('NSE=',num2str(round(NSE,2)))});

