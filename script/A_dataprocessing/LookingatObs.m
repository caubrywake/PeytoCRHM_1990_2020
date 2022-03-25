% load obs
load('D:\PeytoCRHM_1990_2020\data_raw\peyto_met\PeytoMNH_1987_2020.mat')

T = timetable(PeytoMNH19872020_time,PeytoMNH19872020);
TT = retime(T, 'monthly','mean');
timeobs_mth = TT.PeytoMNH19872020_time;
obs_mth = table2array(TT);

T = timetable(PeytoMNH19872020_time,PeytoMNH19872020(:, 6));
TT = retime(T, 'yearly','sum');
timeobs_yr = TT.PeytoMNH19872020_time;
p_yr = table2array(TT);
figure; plot(timeobs_mth, obs_mth(:, 1))
figure; plot(PeytoMNH19872020_time, cumsum(PeytoMNH19872020(:,6)))
figure; bar(timeobs_yr, p_yr)

figure; plot(timeobs_mth, obs_mth(:,4))