%% Model eval 
%% Mass balance

%% This evaulation winter and sumemr mass balance

%% setp up
% close all
% clear all

addpath(genpath('D:\PeytoCRHM_1990_2020\chrm'))
addpath('D:\PeytoCRHM_1990_2020\function')
savedir = 'D:\PeytoCRHM_1990_2020\data_process\modeleval\'
figdir = 'D:\PeytoCRHM_1990_2020\fig\modeleval\'
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\SWE.mat', 'SWE', 'timeCRHM')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev')

% Load validation data
load('D:\FuturePeyto\data\evaluationdata\NB_elevband_WGMS.mat')
load('D:\FuturePeyto\data\evaluationdata\Stake_GSC.mat')
load('D:\FuturePeyto\data\evaluationdata\WB_elevband_WGMS.mat')
load('D:\FuturePeyto\data\evaluationdata\WGMS_summer.mat')
%
%%
fig = figure('units','inches','outerposition',[0 0 8 6]);
yr = [1989:1990, 1993:1995, 2003:2009, 2011:2018];
hru = [1:16, 18:20, 34];

% ice = glacier_h2o;
for i = 1:length(yr)
    subplot(5,4,i)
% GSC
a = find(Stake(:, 1) == yr(i));
sc1 = scatter(Stake(a, 5), Stake(a, 6)/100,10, '^', 'k'); hold on; 
% WGMS
a = find(WGMS_winter(:, 1) == yr(i));
sc1 = scatter(WGMS_winter(a, 2), WGMS_winter(a, 3)/1000,10, '^', 'k'); hold on; 

 % HRU CRHM
 b = find(timeCRHM== datetime(strcat('1-May-', num2str(yr(i)))));
            sc3 = scatter(hruelev(hru), SWE(b, hru)/1000,10, 'd', 'filled', 'b'); hold on
xlim([2100 3000]);
ylim([0 2.2]);
yticks([0.5:0.5:2.2]);
xticks ([2100:250:3050]);

text(2130, 2, num2str(yr(i)))
 a = get(gca,'XTickLabel');  
 set(gca,'XTickLabel',a,'fontsize',10);
  a = get(gca,'YTickLabel');  
  set(gca,'YTickLabel',a,'fontsize',10);
grid on
box on
set(gca, 'FontSize', 10');

end 

% leg = legend ([sc1 sc2 sc3], 'Measured', 'Modelled', 'Modelled (WRF precip offset)',...
%      'Orientation', 'vertical', 'Location', [.5 0.95 0 0] ); 
% ylabel({'Winter Mass Balance';'(m w.e.)'},supAxes)
% xlabel('Elevation (m.a.s.l.)',supAxes)

% 
figname ='WinterMassBalance';
saveas (gcf, strcat( savedir, figname, '.pdf'))
saveas (gcf, strcat(savedir, figname, '.png'))
savefig(gcf, strcat(savedir, figname))
%%

%% find the stanadr deviation and the lapse rate for winter 
% make a graph with some quantifivcation of uncertainty:
% for CRHM: stadnadard devaiation of value at each HRY
% for CRHM, average at each hru
yr = 2003:2018;
hru = [1:9, 16, 18:20, 33];

% find the end of winter swe for individual hru
for  i = 1:length(yr) % select end of winter SWE for each HRU for every year
    b = find(timeCRHM== datetime(strcat('1-May-', num2str(yr(i)))));
for e = 1:length(hru) % for the 30 hr
WinterBalance_CRHM(e, i) = SWE(b, hru(e))/1000;
end 
end 
% Find mean and standard deviation
WinterBalance_CRHM_mean = mean(WinterBalance_CRHM, 2)
WinterBalance_CRHM_std = std(WinterBalance_CRHM, 0, 2)
% Sort by elevation
[elev,id] = sort(hruelev(hru));
WinterBalance_CRHM_mean = WinterBalance_CRHM_mean(id);
WinterBalance_CRHM_std = WinterBalance_CRHM_std (id);
%  error plot
figure;
cr = [0 0 180]/255
errorbar(elev, WinterBalance_CRHM_mean , WinterBalance_CRHM_std , 'vertical','x', 'Color', cr, 'CapSize', 0); hold on
legend ('CRHM HRU winter balance')

% For GSC
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev')
load('D:\PeytoCRHM_1990_2020\data_raw\massbalance\StakePeyto_2005_2020_withdates.mat')
Stake = table2array(StakePeyto20052020withdates(:, 1:8));
Stake_bwdate = table2array(StakePeyto20052020withdates(:, 9));
Stake_bsdate = table2array(StakePeyto20052020withdates(:, 10));

elev_less2900= elev; elev_less2900(elev_less2900>2900) = []; % select elevation below 2900
Sx = Stake(:, 5:6); % select stake winter balance and elevation
WinterBalance_GSC_std  = std(Sx(:, 2))/100;% standard deviation of Stake balance

ck = [150 150 150]./255; % colro for plot
sz = 10; % size of marker
scatter(Sx(:,1), Sx(:, 2)/100, sz, ck, 'filled'); % plot the stake measurement

% find the best fit line to stake winter balance
a = find(isnan(Sx(:, 1)));Sx(a,:) = [];
a = find(isnan(Sx(:, 2)));Sx(a,:) = [];
[pGSC] = polyfit (Sx(:, 1), Sx(:, 2)/100,1) % p(2) is the lapse rate
WinterBalance_bestFit_GSC = polyval(pGSC, elev_less2900);
plot(elev_less2900, WinterBalance_bestFit_GSC, 'k')% add best fit line to plot

% CRHM best fit line
WinterBalanceCRHM_rshp = reshape(WinterBalance_CRHM, 16*14, 1);% all winter balance point in one line
elev_reshape = reshape(repmat(hruelev(hru), 16, 1), 16*14, 1)'; %all elevation in one line
WinterBalance_CRHM_std_allelev= std(WinterBalanceCRHM_rshp)
% Calculate Std without HRU 24 (lots of avalanche activity)
WinterBalanceCRHM_rshp_not34 = WinterBalanceCRHM_rshp; % copy the winter balance matrix
x = find(elev_reshape == hruelev(33)); % fond the leevation corrspsonding to HRU 34
WinterBalanceCRHM_rshp_not34(x) = []; % rmeove these evalues
WinterBalance_CRHM_std_not34 = std(WinterBalanceCRHM_rshp_not34) % not a big difference

% Obtain CRHM best fit line for HRU below 2900
WinterBalance_CRHM_less2900 = WinterBalanceCRHM_rshp ;% copy the winter balance
elev_reshape_less2900 = elev_reshape;
x = find(elev_reshape_less2900>2900); % select the highest elevations
elev_reshape_less2900(x) = [];% remove the highest 2 hru
WinterBalance_CRHM_less2900(x) = [];% remove the highest 2 hru
WinterBalance_CRHM_less2900_std = std(WinterBalance_CRHM_less2900)
[pchrm] = polyfit (elev_reshape_less2900, WinterBalance_CRHM_less2900',1); % pchrm(2) is lapse rate
ychrm = polyval(pchrm, elev_less2900);
% plot(elev_less2900, ychrm, 'b');
legend ('CRHM HRU', 'GSC Stakes', 'GSC best fit line', 'CRHM best fit line')
ylabel ('Winter Mass Balance (m w.e.)');
xlabel ('Elevation (m.a.s.l.')

% % save figure
figname ='WinterBalanceBalance_STD_LapseRate';
saveas (gcf, strcat( savedir, figname, '.pdf'))
saveas (gcf, strcat( savedir, figname, '.png'))
savefig(gcf, strcat( savedir, figname))

%% Summer Mass Balance
% 
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\glacierh2o.mat', 'glacierh2o', 'timeCRHM')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\firn_ice.mat')

% yr = [2003:2009, 2012:2015];
yr = [1993:1995, 2003:2009, 2011:2018];
hru = [1:16, 18:20, 34];
ice = ice+firn;
% how mnay winter mass balance 
Xwinter = Stake( 51:end, 6);
Xwinter = Xwinter(~isnan(Xwinter))
Xsummer = Stake( 51:end, 7);
Xsummer = Xsummer(~isnan(Xsummer))
min(Stake(:, 5))% minimum elevation of stake
max(Stake(:, 5)) % maximum stake elevation
fig = figure('units','inches','outerposition',[0 0 8 6]);

for i = 1:length(yr)
    subplot(5,4,i)
 % GSC
a = find(Stake(:, 1) == yr(i));
scatter(Stake(a, 5), Stake(a, 7)/100, 10, 'k', 'o'); hold on; 

% WGMS
a = find(WGMS_summer(:, 1) == yr(i));
sc1 = scatter(WGMS_summer(a, 2), WGMS_summer(a, 3)/1000,10, '^', 'k'); hold on; 

% CRHM HRU

b = find(timeCRHM== datetime(strcat('1-Apr-', num2str(yr(i)))));
c = find(timeCRHM== datetime(strcat('30-Sep-', num2str(yr(i)))));
MIN1 = ice(c, hru);
MIN2 = ice(b,hru);
summer = (MIN1-MIN2);
summer(summer==0) = nan;
scatter(hruelev(hru), summer/1000,10,  'r', 'd', 'filled' ); hold on

% Formatting
xlim([2100 2850]);
ylim([-7.5 0.7]);
yticks([-5:2.5:0]);
xticks ([2250:250:2750]);
text(2130, 0, num2str(yr(i)))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10);
a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',10);
grid on
box on
set(gca, 'FontSize', 10');

end 

% leg = legend ('Measured', 'Modelled',...
%      'Orientation', 'Horizontal', 'Location', [.5 0.95 0 0]); 
%   ylabel({'Summer Mass Balance';'(m w.e.)'}); %,supAxes)
% xlabel('Elevation (m.a.s.l.)'); %,supAxes)

% tightfig(fig)
figname ='SummerMassBalance';
saveas (gcf, strcat( savedir, figname, '.pdf'))
saveas (gcf, strcat(savedir, figname, '.png'))
savefig(gcf, strcat(savedir, figname))
%% find the stanadr deviation and the lapse rate for summer
yr = 2003:2018;
hru = [1:16, 18:20, 33:34];
% ice = glacierh2o;
% find the end of summer melt  for individual hru
for  i = 1:length(yr) % select end of winter SWE for each HRU for every year
    b = find(timeCRHM== datetime(strcat('1-Apr-', num2str(yr(i)))));
    c = find(timeCRHM== datetime(strcat('15-Sep-', num2str(yr(i)))));
for e = 1:length(hru)
    MIN1 = ice(c, hru(e));
    MIN2 = ice(b,hru(e));
melt = (MIN1-MIN2);
melt(melt==0) = nan;
SumBalance_CRHM(e, i) = melt/1000;
end 
end 
% Find mean and standard deviation
SumBalance_CRHM_mean = nanmean(SumBalance_CRHM, 2)
SumBalance_CRHM_std = nanstd(SumBalance_CRHM, 0, 2)
% Sort by elevation
[elev,id] = sort(hruelev(hru));
SumBalance_CRHM_mean = SumBalance_CRHM_mean(id);
SumBalance_CRHM_std = SumBalance_CRHM_std (id);
%  error plot
figure;
cr = [180 0 0 ]/255
errorbar(elev, SumBalance_CRHM_mean , SumBalance_CRHM_std , 'vertical','x', 'Color', cr, 'CapSize', 0); hold on
legend ('CRHM HRU summer balance')

% For GSC
elev_less2900= elev; elev_less2900(elev_less2900>2900) = []; % select elevation below 2900
Sx = Stake(:, [5,7]); % select stake winter balance and elevation
SumBalance_GSC_std  = nanstd(Sx(:, 2))/100;% standard deviation of Stake balance

ck = [150 150 150]./255; % colro for plot
sz = 10; % size of marker
scatter(Sx(:,1), Sx(:, 2)/100, sz, ck, 'filled'); % plot the stake measurement

% find the best fit line to stake winter balance
a = find(isnan(Sx(:, 1)));Sx(a,:) = [];
a = find(isnan(Sx(:, 2)));Sx(a,:) = [];
[pGSC_sum] = polyfit (Sx(:, 1), Sx(:, 2)/100,1) % p(2) is the lapse rate
SumBalance_bestFit_GSC = polyval(pGSC_sum, elev_less2900);
plot(elev_less2900, SumBalance_bestFit_GSC, 'k')% add best fit line to plot

% CRHM best fit line
SumBalanceCRHM_rshp = reshape(SumBalance_CRHM, 13*20, 1);% all winter balance point in one line
elev_reshape = reshape(repmat(hruelev(hru), 13, 1), 13*20, 1)'; %all elevation in one line
SumBalance_CRHM_std_allelev= nanstd(SumBalanceCRHM_rshp)
% Calculate Std without HRU 24 (lots of avalanche activity)
SumBalanceCRHM_rshp_not34 = SumBalanceCRHM_rshp; % copy the winter balance matrix
x = find(elev_reshape == hruelev(33)); % fond the leevation corrspsonding to HRU 34
SumBalanceCRHM_rshp_not34(x) = []; % rmeove these evalues
SumBalance_CRHM_std_not34 = nanstd(SumBalanceCRHM_rshp_not34) % not a big difference

% Obtain CRHM best fit line for HRU below 2900
SumBalance_CRHM_less2900 = SumBalanceCRHM_rshp ;% copy the winter balance
elev_reshape_less2900 = elev_reshape;
x = find(elev_reshape_less2900>2900); % select the highest elevations
elev_reshape_less2900(x) = [];% remove the highest 2 hru
SumBalance_CRHM_less2900(x) = [];% remove the highest 2 hru
% remove nan
x = find(isnan(SumBalance_CRHM_less2900));
elev_reshape_less2900(x) = [];% 
SumBalance_CRHM_less2900(x) = [];% 
SumBalance_CRHM_less2900_std = nanstd(SumBalance_CRHM_less2900)
[pchrm_sum] = polyfit (elev_reshape_less2900, SumBalance_CRHM_less2900',1); % pchrm(2) is lapse rate

SumBalance_bestFit_CRHM = polyval(pchrm_sum, elev_less2900);
% plot(elev_less2900, SumBalance_bestFit_CRHM, 'r');
legend ('CRHM HRU', 'GSC Stakes', 'GSC best fit line', 'CRHM best fit line', 'location', 'southeast')
ylabel ('Summer Mass Balance (m w.e.)');
xlabel ('Elevation (m.a.s.l.')

% difference in summer melt
SumBalance_BestFit_diff = mean(SumBalance_bestFit_CRHM - SumBalance_bestFit_GSC)
%figure; plot(SumBalance_bestFit_CRHM - SumBalance_bestFit_GSC)
% save figure
figname ='SumBalanceBalance_STD_LapseRate';
saveas (gcf, strcat( savedir, figname, '.pdf'))
saveas (gcf, strcat( savedir, figname, '.png'))
savefig(gcf, strcat( savedir, figname))

%% SWE lapse rate as aspect
% % 
% yr = [2003:2009, 2012:2015];
% fig = figure('units','inches','outerposition',[0 0 9 7]);
% 
% subplot(2,2,1)
% hru = [1:10];
% for i = 1:11
%     a = find(Stake(:, 1) == yr(i));
%     stakeyr = [Stake(a, 5), Stake(a, 6)*10] ;
%     [StyrSrt, idsort] = sort(stakeyr(:, 1));
%     WBsrt = stakeyr(idsort, 2);
%     plot(StyrSrt, WBsrt, '-ok'); hold on
%     
%     b = find(timeCRHM== datetime(strcat('1-May-', num2str(yr(i)))));
%   crhmyr = [hruelev(hru), SWECUR(b, hru)'] ;
%     [HRUyrSrt, idsortcrhm] = sort(crhmyr(:, 1));
%     SWEsrt = crhmyr(idsortcrhm, 2);
%     plot(HRUyrSrt, SWEsrt, '-xr'); hold on
%         ylim ([0 1800])
% end 
% title ('Flowline 1, HRU 1-10 (North facing)')
% ylabel ('Winter Mass Balance (mm w.e.)')
% xlabel('Elevation (m.a.s.l.)')
% legend ('NRCan Stake', 'CRHM HRU', 'Location', 'NorthWest')
% subplot(2,2,2)
% hru = [10:17];
% for i = 1:11
%     a = find(Stake(:, 1) == yr(i))
%     stakeyr = [Stake(a, 5), Stake(a, 6)*10] 
%     [StyrSrt, idsort] = sort(stakeyr(:, 1))
%     WBsrt = stakeyr(idsort, 2)
%     plot(StyrSrt, WBsrt, '-ok'); hold on
%     
%     b = find(timeCRHM== datetime(strcat('1-May-', num2str(yr(i)))));
%   crhmyr = [hruelev(hru), SWECUR(b, hru)'] 
%     [HRUyrSrt, idsortcrhm] = sort(crhmyr(:, 1))
%     SWEsrt = crhmyr(idsortcrhm, 2)
%     plot(HRUyrSrt, SWEsrt, '-xr'); hold on
%         ylim ([0 1800])
% end 
% title ('Flowline 2, HRU 11-17 (East facing)')
% ylabel ('Winter Mass Balance (mm w.e.)')
% xlabel('Elevation (m.a.s.l.)')
% 
% subplot(2,2,3)
% hru = [18:21];
% for i = 1:11
%     a = find(Stake(:, 1) == yr(i))
%     stakeyr = [Stake(a, 5), Stake(a, 6)*10] 
%     [StyrSrt, idsort] = sort(stakeyr(:, 1))
%     WBsrt = stakeyr(idsort, 2)
%     plot(StyrSrt, WBsrt, '-ok'); hold on
%     
%     b = find(timeCRHM== datetime(strcat('1-May-', num2str(yr(i)))));
%   crhmyr = [hruelev(hru), SWECUR(b, hru)'] 
%     [HRUyrSrt, idsortcrhm] = sort(crhmyr(:, 1))
%     SWEsrt = crhmyr(idsortcrhm, 2)
%     plot(HRUyrSrt, SWEsrt, '-xr'); hold on
%         ylim ([0 1800])
% 
% end 
% title ('Flowline 3, HRU 18-22 (South facing)')
% ylabel ('Winter Mass Balance (mm w.e.)')
% xlabel('Elevation (m.a.s.l.)')
% 
% subplot(2,2,4)
% hru = [1:21];
% for i = 1:11
%     a = find(Stake(:, 1) == yr(i))
%     stakeyr = [Stake(a, 5), Stake(a, 6)*10] 
%     [StyrSrt, idsort] = sort(stakeyr(:, 1))
%     WBsrt = stakeyr(idsort, 2)
%     plot(StyrSrt, WBsrt, '-ok'); hold on
%     
%     b = find(timeCRHM== datetime(strcat('1-May-', num2str(yr(i)))));
%   crhmyr = [hruelev(hru), SWECUR(b, hru)'] 
%     [HRUyrSrt, idsortcrhm] = sort(crhmyr(:, 1))
%     SWEsrt = crhmyr(idsortcrhm, 2)
%     plot(HRUyrSrt, SWEsrt, '-xr'); hold on
%     ylim ([0 1800])
% end 
% title ('All Glacier HRU')
% ylabel ('Winter Mass Balance (mm w.e.)')
% xlabel('Elevation (m.a.s.l.)')
% 
% tightfig(fig)
% 
% figname ='WinterMassBalanceGradient_PerFlowline_1OBS';
% saveas (gcf, strcat( 'D:\FuturePeyto\fig\B1b\', figname, '.pdf'))
% saveas (gcf, strcat('D:\FuturePeyto\fig\B1b\', figname, '.png'))
% savefig(gcf, strcat('D:\FuturePeyto\fig\B1b\', figname))
% 
