%% B2 SWE with snow piallow at wildcat and Bow summit surveys
% wildcat snow pillow was established in 2012, at an elevation of 2122
savedir = 'D:\PeytoCRHM_1990_2020\data_process\modeleval\'
figdir = 'D:\PeytoCRHM_1990_2020\fig\modeleval\'

%% compare baisn SWE with snow pillow at wildcat
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat', 'hruelev')
load('D:\PeytoCRHM_1990_2020\data_raw\snowsurvey_bowsummit\BowSummit_AlbertaPark_SWE_1980_2020.mat')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\SWE.mat', 'SWE', 'timeCRHM')
SWE_BSsurvey = table2array(BowSummitAlbertaParkSWE19802020(:, 2));
SWE_BSsurvey_time = table2array(BowSummitAlbertaParkSWE19802020(:, 1));
T = timetable(SWE_BSsurvey_time,SWE_BSsurvey);
TT = retime(T, 'monthly', 'mean');
SWE_BSsurvey_time = TT.SWE_BSsurvey_time;
SWE_BSsurvey = table2array(TT);
T = timetable(timeCRHM,SWE);
TT = retime(T, 'monthly', 'mean');
SWEm = table2array(TT);
timeCRHMm= TT.timeCRHM;

figure
plot(SWE_BSsurvey_time, SWE_BSsurvey, '-xk'); hold on
plot(timeCRHMm, SWEm(:, 8), '-xr')
plot(timeCRHMm, SWEm(:, 9), '-xr')
plot(timeCRHMm, SWEm(:, 10), '-xr')
ylim([0 1500])

%% wildcat: 2122 m (6962 ft)     51.696° N, 116.629° W
load('D:\PeytoCRHM_1990_2020\data_raw\snowpillow_wildcat\DataSetExport-SW.Daily@2A32P-20210305233309.mat')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\SWE.mat', 'SWE', 'timeCRHM')

 % https://aqrt.nrs.gov.bc.ca/Data/Location/Summary/Location/2A32P/Interval/Latest
SWE_wildcat = table2array(DataSetExportSW(:, 2));
SWE_wildcat_time = datetime(table2array(DataSetExportSW(:, 1)));

fig = figure('units','inches','outerposition',[0 0 7 4]);
fig.Renderer='Painters';
p1 = plot(timeCRHM, SWE(:, 8), 'Color', [57 206 215]./255, 'linewidth' , 0.6);hold on
p2 = plot(timeCRHM, SWE(:, 17), 'Color',  [25 199 95]./255, 'linewidth' , 0.6); 
p3 = plot(timeCRHM, SWE(:, 7), 'Color',  [25 112 221]./255, 'linewidth' , 0.6); 
p4 = plot(SWE_wildcat_time, SWE_wildcat, 'k', 'linewidth' , 1); hold on
xlim([datetime('01-Oct-2014') datetime('01-Oct-2020')])
legend ([p4(1) p1(1) p2(1) p3(1)],  'Wildcat SP, 2122m','CRHM, 2211m', 'CRHM, 2251m',  'CRHM, 2252m','location', 'Northoutside', 'orientation', 'horizontal')
ylabel ('SWE (mm w.e.)')
xticks ([datetime('01-Jan-2015'):calmonths(6):datetime('01-Oct-2020')])
xtickformat ('dd-MM-yyyy')
xtickangle (30)
figname ='SWE_2014_2020_SnowPillow';
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'), '.pdf'))
saveas (gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM'),  '.png'))
savefig(gcf, strcat(figdir, figname, datestr(now, 'yyyymmddHHMM')))


t1 = find(timeCRHM == SWE_wildcat_time(1));
 t = timeCRHM(t1):days(1):datetime('28-Sep-2020 16:00:00');
T = timetable(SWE_wildcat_time, SWE_wildcat);
TT = retime(T, t, 'linear');
t_wildcat = TT.SWE_wildcat_time;
wildcat = table2array(TT);
T = timetable(timeCRHM, SWE);
TT = retime(T, t, 'linear');
t_chrm = TT.timeCRHM;
crhm = table2array(TT);
figure; plot(t_wildcat, wildcat)
hold on; plot(t_chrm, crhm(:, 7))
corrcoef(wildcat, crhm(:, 8))
corrcoef(wildcat, crhm(:, 17))
corrcoef(wildcat, crhm(:, 7))

