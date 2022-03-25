%% Export SR50 as obs file
close all
clear all
load('D:\SaltDilution\F_buildingratingcurve\output\Q_SR50_daily_20132020.mat')
t = Q_SR50_daily_time+days(1);
t = datevec(t);
t(:, 4)=1;

sr50_obs= [t(1:end, 1:5) Q_SR50_daily(1:end, :)*3600]; % compiled time and data together, and save it in a matlab format

% create the obs text file
headerlines = {'Obs file, CUR bc';
              'alb 1 ()';
              '#####'}
fp = strcat('D:\PeytoCRHM_1990_2020\data_raw\sr50_water\sr50_water.obs');
fid = fopen(fp, 'wt');
for l = 1:numel(headerlines)
   fprintf(fid, '%s\n', headerlines{l});
end
fclose(fid);
dlmwrite(fp, sr50_obs , '-append', 'delimiter', '\t'); 

