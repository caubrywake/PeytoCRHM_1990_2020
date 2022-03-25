%% B3 
%% evaluation of the SR50
%% for three location: 
% lower ice, middle ice, upper ice
% elevation 

% lower ice: 2173-2183m - hru 7-8
% middle ice = 2454-2461 - hru 6
% upper ice = 2709 - hru3
% load('PeytoMiddleICE_SR50_hrly_17Sept2010_14Sept2013.mat')
% load('PeytoUpperICE_SR50_hrly_5Aug2011_5Sept2013.mat')
close all
clear all

figdir = 'D:\4_PeytoCRHM_1990_2020\fig\modeleval\'

load('D:\4_PeytoCRHM_1990_2020\data_raw\sr50_ice\SR50_LowerIce_2010_2020.mat')

load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\glacierh2o.mat')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\firn_ice.mat')
%% create time array
% upper_t = datetime([PeytoUpper(:, 1:5) zeros(length(PeytoUpper),1)]);
% upper = PeytoUpper(:, 6);
% middle_t = datetime([PeytoMiddle(:, 1:5) zeros(length(PeytoMiddle),1)]);
% middle = PeytoMiddle(:, 6);
lower = SR50_LowerIce_2010_2020;
lower_t = SR50_LowerIce_2010_2020_time;
lower(70385:71349)=[];
lower_t(70385:71349)=[];
% plot(upper_t, upper);hold on
% plot(middle_t, middle);hold on

% change lower to mm w.e. 
lower = -lower; % now its melting
lower = lower.*1000; % in mm
lower = lower.*900./1000;
figure
plot(lower_t, lower(:,1));

%% Compare modelled SWE and on ice sr50
% % select all the time step where ther can be comparison
% sr50_validtime = [datetime('04-May-2011 01:00') datetime('05-Oct-2011 16:00'); ...
%      datetime('19-Apr-2012 01:00') datetime('8-Aug-2012 16:00'); ...
%      datetime('25-Aug-2012 21:00') datetime('15-Oct-2012 7:00'); ... 
%      datetime('2-May-2013 12:00') datetime('11-Jul-2013 13:00');...
%      datetime('14-Aug-2013 21:00') datetime('20-Sep-2013 23:00'); ...
%      datetime('6-May-2014 10:00') datetime('27-Sep-2014 01:00'); ...
%      datetime('3-May-2015 10:00') datetime('12-Jun-2015 11:00'); ...
%      datetime('22-Jul-2015 7:00') datetime('31-Aug-2015 7:00'); ...
%      datetime('28-Mar-2016 16:00') datetime('10-Oct-2016 17:00'); ...
%      datetime('2-May-2017 16:00') datetime('27-Jun-2017 10:00'); ...
%      datetime('9-Aug-2017 23:00') datetime('18-Sep-2017 9:00'); ...
%      datetime('21-Apr-2018 14:00') datetime('15-May-2018 10:00'); ...
%      datetime('01-Aug-2018 15:00') datetime('29-Sep-2018 21:00'); ...
%      datetime('04-May-2019 15:00') datetime('27-Aug-2019 05:00')];
%  
 %% part 2
 sr50_validtime = ...
     [datetime('5-May-2011 01:00') datetime('05-Oct-2011 16:00'); ...
     datetime('8-Jun-2012 013:00') datetime('8-Aug-2012 16:00'); ...
     datetime('25-Aug-2012 21:00') datetime('15-Oct-2012 7:00'); ... 
     datetime('2-May-2013 12:00') datetime('11-Jul-2013 13:00');...
     datetime('14-Aug-2013 21:00') datetime('20-Sep-2013 23:00'); ...
     datetime('23-Jun-2014 10:00') datetime('27-Sep-2014 01:00'); ...
     datetime('19-May-2015 10:00') datetime('12-Jun-2015 11:00'); ...
     datetime('22-Jul-2015 7:00') datetime('31-Aug-2015 7:00'); ...
     datetime('28-Mar-2016 16:00') datetime('10-Oct-2016 17:00'); ...
     datetime('25-May-2017 16:00') datetime('27-Jun-2017 10:00'); ...
     datetime('9-Aug-2017 23:00') datetime('18-Sep-2017 9:00'); ...
     datetime('21-Apr-2018 14:00') datetime('15-May-2018 10:00'); ...
     datetime('01-Aug-2018 15:00') datetime('29-Sep-2018 21:00'); ...
     datetime('5-Jun-2019 15:00') datetime('27-Aug-2019 05:00')];

  for i = 1:length(sr50_validtime)
sr50_validtime_CRHM(i, 1) = find(timeCRHM== sr50_validtime(i, 1));
sr50_validtime_CRHM(i, 2) = find(timeCRHM== sr50_validtime(i, 2));
sr50_validtime_meas(i, 1) = find(lower_t== sr50_validtime(i,1));
sr50_validtime_meas(i, 2) = find(lower_t== sr50_validtime(i,2));
  end 
  
  %% Figure of each event
close all
clear R RMSE NSE MBE melttot

fig = figure('units','inches','outerposition',[0 0 8 10]);
fig.Renderer='Painters';
 ice2 = glacierh2o(:, 7);
  ice3 = glacierh2o(:, 8);
  lab = {'(a)';'(b)';'(c)';'(d)';'(e)';'(f)';'(g)';'(h)';'(i)';'(j)';'(k)';'(l)';'(m)';'(n)'}
 for i = 1:length(sr50_validtime_meas)
ice_mod = ice2(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2)) - ice2(sr50_validtime_CRHM(i, 1)); 
ice_meas = lower(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) - lower(sr50_validtime_meas(i, 1)); 
t = timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2));

T = timetable(t, ice_meas, ice_mod);
TT = retime(T, 'daily','mean');
ice_mod = table2array(TT(:,2));
ice_meas = table2array(TT(:,1));
t = TT.t;

% statistics
r = corrcoef(ice_meas, ice_mod);
R(i) = r(2);
RMSE(i) = sqrt(mean((ice_meas - ice_mod).^2));
NSE(i) =  nashsutcliffe([datenum(t),ice_meas],  [datenum(t),ice_mod]);
MBE(i) =  (ice_mod(end)-ice_meas(end))/ice_meas(end);
totmelt (i, 1) = ice_meas(end);
totmelt(i, 2) = ice_mod(end);
nb_day(i) = numel(t);
subplot(4,4,i)
plot(t, ice_mod/1000, 'r'); hold on
% plot(timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2)), ice_mod2, 'm'); 
hold on
plot (t, ice_meas/1000 , 'k')
if min(ice_meas)<=min(ice_mod)
ylim([min(ice_meas)/1000 max(ice_mod)/1000]); 
else
   ylim([min(ice_mod)/1000 max(ice_mod)/1000])
end 
xlim ([timeCRHM(sr50_validtime_CRHM(i, 1)) timeCRHM(sr50_validtime_CRHM(i, 2))]);
if i == 1
lg = legend ('modelled (2252 m)', 'measured (2173-2183m)', 'orientation', 'horizontal', 'location', 'northeast')
end 
TextLocation(lab(i),'Location','NorthWest');

% text
str1 = strcat('NSE = ', num2str(round(NSE(i),2))); 
str2 = strcat('MB = ', num2str(round(MBE(i),2))); 
str3 = {str1;str2};
TextLocation(str3,'Location','SouthWest');
xtickformat('dd-MMM-yy')
xtickangle(25)
 end
%  pos = get(lg,'Position');
% set(lg,'Position',[pos(1)+0.095 pos(2)+0.075 pos(3) pos(4)]);

 [ax2, h2] = suplabel('Surface Melt (m w.e.)', 'y');
 ax2.Position = ax2.Position + [+0.03 0 0 0];

 tightfig(fig)
 figname ='SR50_eval_panel';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))

%% statistic for the whole period together?
clear R RMSE NSE MBE melttot
 for i = 1:length(sr50_validtime_meas)
ice_mod = ice2(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2)) - ice2(sr50_validtime_CRHM(i, 1)); 
ice_meas = lower(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) - lower(sr50_validtime_meas(i, 1)); 
t = timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2));

T = timetable(t, ice_meas, ice_mod);
TT = retime(T, 'daily','mean');
ice_mod = table2array(TT(:,2));
ice_meas = table2array(TT(:,1));
t = TT.t;
 
if i ==1
MEAS = ice_meas;
MOD = ice_mod;
tt = t;
else 
MEAS = [MEAS; ice_meas];
MOD = [MOD; ice_mod];
tt = [tt;t];
end 
 end 
% statistics
r = corrcoef(MEAS, MOD);
R_all= r(2);
RMSE_all = sqrt(mean((MEAS - MOD).^2))/1000;
NSE_all =  nashsutcliffe([datenum(tt),MEAS],  [datenum(tt),MOD]);
MBE_all =  (MOD(end)-MEAS(end))/MEAS(end);
%% statistics
% cumulative melt difference


%% All in one -continuous
fig = figure('units','inches','outerposition',[0 0 7 3]);
fig.Renderer='Painters';

ice2 = glacierh2o(:, 7);
ice2 = ice2 - ice2(sr50_validtime_CRHM(1, 1));
ice2 = ice2(sr50_validtime_CRHM(1, 1):sr50_validtime_CRHM(14, 2))/1000;
tall  = timeCRHM(sr50_validtime_CRHM(1, 1):sr50_validtime_CRHM(14, 2));
T = timetable(tall, ice2);
TT = retime(T, 'daily','mean');
ice2 = table2array(TT);
t2 = TT.tall;
plot(t2, ice2,'r')
hold on
for i = 1:length(sr50_validtime_meas)
ice_meas = lower(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) - lower(sr50_validtime_meas(i, 1)); 
t = timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2));
T = timetable(t, ice_meas);
TT = retime(T, 'daily','mean');
ice_meas = table2array(TT)/1000;
t = TT.t;

t1 = find(t2 == t(1))
ice_meas_delta = ice_meas +(ice2(t1)- ice_meas(1)); 

plot(t, ice_meas_delta, 'k', 'linewidth', 2)

end 
ylabel ('Surface Melt (m w.e.)')
legend ('Modelled Surface Melt','Measured Surface Melt') 
ylim ([-50 0.1])
xlim([timeCRHM(sr50_validtime_CRHM(1, 1))-days(30) datetime('01-Oct-2019')])
 tightfig(fig)
 figname ='SR50_eval_continous';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))

%%
figure
ice = glacierh2o(:, 7);
ice = ice- ice(sr50_validtime_CRHM(1, 1));
plot(timeCRHM(sr50_validtime_CRHM(1, 1):sr50_validtime_CRHM(14, 2)), ice(sr50_validtime_CRHM(1, 1):sr50_validtime_CRHM(14, 2)),'k')
hold on
for i = 1:length(sr50_validtime_meas)
ice_meas = lower(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) ;
%plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), ice_meas)

ice_meas_delta = ice_meas +(ice(sr50_validtime_CRHM(i, 1))- lower(sr50_validtime_meas(i, 1))); 
plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), ice_meas_delta, 'r')
end 


