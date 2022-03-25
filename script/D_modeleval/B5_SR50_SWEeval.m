%% B4 - SWE eval using SR50
% redo but with snow depth simulation

% load('PeytoUpperICE_SR50_hrly_5Aug2011_5Sept2013.mat')
close all
clear all

figdir = 'D:\4_PeytoCRHM_1990_2020\fig\modeleval\'
load('D:\4_PeytoCRHM_1990_2020\data_raw\sr50_ice\SR50_LowerIce_2010_2020.mat')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\SWE.mat')

 delimiterIn = '\t';
 [D, delimiterOut, headerlinesOut]= importdata('D:\4_PeytoCRHM_1990_2020\chrm\output\cRHMsim_snowdepth.obs', delimiterIn);
 s = D.textdata;
 ss = s(headerlinesOut+1:end, 1);
timeCRHM = datetime(ss, 'InputFormat', 'yyyy MM dd HH mm');
snowdepth = D.data;

lower = SR50_LowerIce_2010_2020;
lower_t = SR50_LowerIce_2010_2020_time;
lower(70385:71349)=[];
lower_t(70385:71349)=[];
% plot(upper_t, upper);hold on
% plot(middle_t, middle);hold on
 a = find(lower_t == '14-Feb-2018 9:00');
 b = find(lower_t == '12-Feb-2018 14:00');
 lower(a:end) = lower(a:end)+(lower(b)-lower(a));

% still in m 
lower = -lower; % now its melting
figure
plot(lower_t, lower(:,1));

% fill in the gap of sr50
T = timetable(lower_t, lower);
TT = retime(T, lower_t, 'linear');
lower_fill = table2array(TT);
figure
plot(lower_t, lower_fill(:,1)); hold on
 plot(lower_t, lower(:,1));
 
%% make continu where there are reset
figure
plot(lower_fill(:,1));

clear SR50_LowerIce_2010_2020_time SR50_LowerIce_2010_2020

%% Select valid winters
 sr50_validtime = ...
     [datetime('10-Oct-2010 21:00') datetime('01-Jul-2011 00:00'); ...
      datetime('10-Oct-2011 13:00') datetime('01-Jul-2012 00:00'); ...
      datetime('16-Oct-2012 00:00') datetime('01-Jul-2013 00:00'); ... 
      datetime('25-Oct-2013 01:00') datetime('01-Jul-2014 00:00');...
      datetime('01-Oct-2014 00:00') datetime('01-Jul-2015 10:00'); ...
      datetime('26-Oct-2015 00:00') datetime('01-Jul-2016 00:00');...
      datetime('04-Oct-2016 00:00') datetime('01-Jul-2017 00:00');...
      datetime('06-Oct-2017 00:00') datetime('10-May-2018 00:00'); ...
      datetime('12-Sep-2018 00:00') datetime('01-Jul-2019 00:00')];
 

  for i = 1:length(sr50_validtime)
sr50_validtime_CRHM(i, 1) = find(timeCRHM== sr50_validtime(i, 1));
sr50_validtime_CRHM(i, 2) = find(timeCRHM== sr50_validtime(i, 2));
sr50_validtime_meas(i, 1) = find(lower_t== sr50_validtime(i,1));
sr50_validtime_meas(i, 2) = find(lower_t== sr50_validtime(i,2));
  end 
%% plot each year
close all

fig = figure('units','inches','outerposition',[0 0 8 10]);
fig.Renderer='Painters';
 ice4= snowdepth (:, 3);
 
  lab = {'(a) 2010';'(b) 2011';'(c) 2012';'(d) 2013';'(e) 2014';'(f) 2015';'(g) 2016';'(h) 2017';'(i) 2018'}
  clear NSE RMSE R
 for i = 1:length(sr50_validtime_meas)
ice_mod3 = ice4(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2)) - ice4(sr50_validtime_CRHM(i, 1)); 

ice_meas = lower_fill(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) - lower_fill(sr50_validtime_meas(i, 1)); 
subplot(3,3,i)
plot(timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2)), ice_mod3, 'b'); hold on

% plot(timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2)), ice_mod2, 'm'); 
hold on
plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), smooth(ice_meas,6) , 'k')
% if min(ice_meas)<=min(ice_mod)
 ylim([0 2.5])
% else
%    ylim([min(ice_mod)/1000 max(ice_mod)/1000])
% end 
 xlim ([timeCRHM(sr50_validtime_CRHM(i, 1)) timeCRHM(sr50_validtime_CRHM(i, 2))]);
 if i == 1
lg = legend ('CRHM, 2176m', 'Measured, 2173-2183m', 'orientation', 'horizontal', 'location', 'northeast')
 pos = lg.Position;
 lg.Position = [ 0.1323    0.9326    0.3942    0.0305];
 end 
xtickformat('dd-MMM')
xtickangle(25)
t1 = datevec(lower_t(sr50_validtime_meas(i, 1)));
yr = t1(:,1);
tickstr = strcat('01-Oct-', num2str(yr));
ticklab = datetime(tickstr)
xticks([ticklab:calmonths(2):ticklab+calmonths(9)])


r = corrcoef(ice_mod3, ice_meas);
R2(i,1) = r(2)^2;
RMSE(i,1) = sqrt(mean((ice_meas - ice_mod3).^2));
NSE(i) =  nashsutcliffe([datenum(timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2))),ice_meas],  [datenum(timeCRHM(sr50_validtime_CRHM(i, 1):sr50_validtime_CRHM(i, 2))),ice_mod3])
MBE(i) =  mean(ice_mod3-ice_meas);
text( timeCRHM(sr50_validtime_CRHM(i, 1))+ days(3), 2.3 , lab{i})
text( timeCRHM(sr50_validtime_CRHM(i, 1))+ days(3), 2.0 , strcat('NSE = ',num2str(round(NSE(i),2))), 'FontSize',8)
text( timeCRHM(sr50_validtime_CRHM(i, 1))+ days(3), 1.7 , strcat('R^2 = ',num2str(round(R2(i),2))), 'FontSize',8)
text( timeCRHM(sr50_validtime_CRHM(i, 1))+ days(3), 1.4, strcat('RMSE = ',num2str(round(RMSE(i),2))), 'FontSize',8)
 end
 
%  pos = get(lg,'Position');
% set(lg,'Position',[pos(1)+0.095 pos(2)+0.075 pos(3) pos(4)]);

 [ax2, h2] = suplabel('Snow depth (m)', 'y');
 ax2.Position = ax2.Position + [+0.03 0 0 0];

 figname ='SR50_SWE_eval_panel';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))
%% elev of hru 7,8,9 and compared to weather station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% continuous
fig = figure('units','inches','outerposition',[0 0 7 4]);
fig.Renderer='Painters';

 ice4= snowdepth (:, 3);
 plot(timeCRHM, smooth(ice4,12), 'b'); hold on


 for i = 1:length(sr50_validtime_meas)

ice_meas = lower_fill(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) - lower_fill(sr50_validtime_meas(i, 1)); 
plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), smooth(ice_meas,12) , 'k')
 end
  xlim([datetime('01-Oct-2010') datetime('01-Jul-2019')]);
  ylim([0 2.5])
  
  %  pos = get(lg,'Position');
% set(lg,'Position',[pos(1)+0.095 pos(2)+0.075 pos(3) pos(4)]);
 lg = legend ('CRHM, 2176m', 'Measured, 2173-2183m', 'orientation', 'horizontal', 'location', 'northeast')
 pos = lg.Position;
 lg.Position = [ 0.4265    0.7982    0.4670    0.0988];


% else
%    ylim([min(ice_mod)/1000 max(ice_mod)/1000])
xtickformat('dd-MMM-yyy')
xtickangle(25)
ylabel('Snow depth (m)');

 figname ='SR50_SWE_eval';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))
%% elev of hru 7,8,9 and compared to weather station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




plot(timeCRHM(sr50_validtime_CRHM(1, 1):sr50_validtime_CRHM(9, 2)), ice4(sr50_validtime_CRHM(1, 1):sr50_validtime_CRHM(9, 2)),'m')
plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), smooth(ice_meas,6) , 'k')

hold on
for i = 1:length(sr50_validtime_meas)
ice_meas = lower_fill(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)) ;
%plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), ice_meas)
ice_meas_delta = ice_meas +(ice2(sr50_validtime_CRHM(i, 1))- lower(sr50_validtime_meas(i, 1))); 
ice_meas_delta(ice_meas_delta<=0)=0;
plot (lower_t(sr50_validtime_meas(i, 1):sr50_validtime_meas(i, 2)), smooth(ice_meas_delta,12), 'k', 'linewidth', 0.8)
end 

ylabel ('Snow Depth (m)')
legend ('Modelled 2252m', 'Modelled 2211m', 'Modelled 2176m', 'Meas 2173-2183m') 
ylim ([0 3])
xlim([timeCRHM(sr50_validtime_CRHM(1, 1))-days(30) datetime('01-Oct-2019')])
 tightfig(fig)
 figname ='SR50_eval_continous';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))

5% for each seaosn, caluclate rmse, r2, mb
