%% make sure there ar enot sesaosnal or montlhy trends
clear all
% compile seasonal anomaly
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWE', 'hru_rain', 'hru_snow', 'hru_t', 'time')
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'basinflow', 'basingw')

savedir = 'D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\'
figdir = 'D:\4_PeytoCRHM_1990_2020\fig\analysis\'
 
%% Compile annual values for annual and seasonal hydrometeorological variables
% for fall, only to 2019
% for for else, to 2020

snow = hru_snow*1000;
rain = hru_rain*1000;
precip = snow+rain;
Ta=hru_t;
timeCRHM = time;
hru_area=hruarea;
hru_elev = hruelev;
yr = 1989:2020;
hru = 1:37;

varname = {'Tw' 'Tsp' 'Tsu' 'Tf' ...
           'Sw' 'Ssp' 'Ssu' 'Sf' ...
           'Rw' 'Rsp' 'Rsu' 'Rf'  ...
           'Pw' 'Psp'  'Psu' 'Pf' ...
           'Ta2' 'Sa' 'Ra' 'Pa' ...
          'Flowsp' 'Flowsu' 'Flowf' 'Flowa'};

%% Montlhy averages
[Tm,~, timem] = MonthlyHruWeigthedBasinSum(timeCRHM,Ta, 1:37, hru_area); 
 [Pm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,precip, 1:37, hru_area); 
[Sm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,snow, 1:37, hru_area); 
[Rm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,rain, 1:37, hru_area); 
[Bm, ~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,basinflow/3600, 1, hru_area);


% Seasonal components
for i = 1:length(yr)
% anuual
   t1 = strcat ('01-Oct-', num2str(yr(i)-1), {' '}, '00:00');
   t2 = strcat ('01-Sep-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Ta2(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Sa(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Ra(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Pa(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flowa (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% Fall (Oct-Dec)
   t1 = strcat ('01-Sep-', num2str(yr(i)-1), {' '}, '00:00');
   t2 = strcat ('01-Nov-', num2str(yr(i)-1), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tf(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Sf(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rf(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Pf(i, 1) = nansum(Pm (a:b)); %m w.e per hru
Flowf (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% Winter (Dec-marc)
    t1 = strcat ('01-Dec-', num2str(yr(i)-1), {' '}, '00:00');
   t2 = strcat ('01-Feb-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem ==datetime(t2));
Tw(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Sw(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rw(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Pw(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flow (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% Spring : Apr-July
   t1 = strcat ('01-Mar-', num2str(yr(i)), {' '}, '00:00');
   t2 = strcat ('01-May-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsp(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Ssp(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rsp(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Psp(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flowsp (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% summer
   t1 = strcat ('01-Jun-', num2str(yr(i)), {' '}, '00:00');
   t2 = strcat ('01-Aug-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsu(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Ssu(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rsu(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Psu(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flowsu (i, 1) = nansum(Bm (a:b)); %m w.e per hru
end 

VarMod= [Tw Tsp Tsu Tf ...
           Sw Ssp Ssu Sf ...
           Rw Rsp Rsu Rf  ...
           Pw Psp  Psu Pf ...
           Ta2 Sa Ra Pa ...
          Flowsp Flowsu Flowf Flowa]
sz = size(VarMod);

% VarMod (10, :) = nan; 

for i = 1:sz(2)
[H(i), p(i)] = Mann_Kendall (VarMod (:, i), 0.05);
[~, ~, h(i), sig(i), ~, ~, sigma(i), sen(i),~, ~, ~, ~, ~, ~, ~, ~] = ktaub([yr', VarMod(:, i)], 0.05, 1)
end
sign_var = find(H);
varname_sign = varname (find(H))
sign_slope = sen(sign_var);
sig_sig =sig(sign_var);

close all
fig = figure('units','inches','outerposition',[0 0 7 5]);

subplot(2,2, 1);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(1))], 0.05, 1);
ylabel ('Summer Ta ({^\circ}C)')
legend ('data', 'Sens slope', 'Slope confidence interval (95%)')

subplot(2,2, 2);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(2))], 0.05, 1);
ylabel ({'Winter Snowfall (mm)'})

subplot(2,2,3);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(3))], 0.05, 1);
ylabel ({'Summer Rainfall (mm)'})

subplot(2,2, 4);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(6))], 0.05, 1);
ylabel ({'Annual Rainfall (mm)'})



filename = 'SignificantTrends_1989_2020'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))

%% Low level values
[Tm,~, timem] = MonthlyHruWeigthedBasinSum(timeCRHM,Ta, [7:10,17, 21, 25,26, 29], hru_area); 
[Pm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,precip,[7:10,17, 21, 25,26, 29], hru_area); 
[Sm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,snow, [7:10,17, 21, 25,26, 29], hru_area); 
[Rm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,rain, [7:10,17, 21, 25,26, 29], hru_area); 
[Bm, ~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,basinflow/3600, 1, hru_area);

% Seasonal components
for i = 1:length(yr)
   t1 = strcat ('01-Dec-', num2str(yr(i)-1));
   t2 = strcat ('01-Mar-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem ==datetime(t2));
Tw(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sw(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rw(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pw (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flow (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Apr-', num2str(yr(i)));
   t2 = strcat ('01-Jul-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsp(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Ssp(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rsp(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Psp (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowsp (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Jul-', num2str(yr(i)));
   t2 = strcat ('01-Oct-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsu(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Ssu(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rsu(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Psu (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowsu (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Oct-', num2str(yr(i)));
   t2 = strcat ('01-Dec-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tf(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sf(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rf(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pf (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowf (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Oct-', num2str(yr(i)-1));
   t2 = strcat ('01-Oct-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Ta2(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sa(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Ra(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pa (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowa (i, 1) = sum(Bm (a:b)); %m w.e per hru
end 

VarMod = [Tw Tsp Tsu Tf ...
    Sw Ssp Ssu Sf...
    Rw Rsp Rsu Rf ...
    Pw Psp Psu Pf...
    Ta2 Sa Ra Pa ...
    Flowsp Flowsu Flowf Flowa]
sz = size(VarMod);

for i = 1:sz(2)
[H(i), p(i)] = Mann_Kendall (VarMod (:, i), 0.05);
[~, ~, h(i), sig(i), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ktaub([VarMod(:, i), yr'], 0.05, 1)
end 
plot(H)
%% High values
clear VarMod Tm Pm Sm Rm Bm
[Tm,~, timem] = MonthlyHruWeigthedBasinSum(timeCRHM,Ta, [1,2,3,10,11,17,18], hru_area); 
[Pm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,precip, [1,2,3,10,11,17,18], hru_area); 
[Sm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,snow, [1,2,3,10,11,17,18], hru_area); 
[Rm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,rain, [1,2,3,10,11,17,18], hru_area); 
[Bm, ~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,basinflow/3600, 1, hru_area);

% Seasonal components
for i = 1:length(yr)
   t1 = strcat ('01-Dec-', num2str(yr(i)-1));
   t2 = strcat ('01-Mar-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem ==datetime(t2));
Tw(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sw(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rw(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pw (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flow (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Apr-', num2str(yr(i)));
   t2 = strcat ('01-Jul-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsp(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Ssp(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rsp(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Psp (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowsp (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Jul-', num2str(yr(i)));
   t2 = strcat ('01-Oct-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsu(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Ssu(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rsu(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Psu (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowsu (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Oct-', num2str(yr(i)));
   t2 = strcat ('01-Dec-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tf(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sf(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rf(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pf (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowf (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Oct-', num2str(yr(i)-1));
   t2 = strcat ('01-Oct-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Ta2(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sa(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Ra(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pa (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowa (i, 1) = sum(Bm (a:b)); %m w.e per hru
end 

VarMod = [Tw Tsp Tsu Tf...
    Sw Ssp Ssu Sf ...
    Rw Rsp Rsu Rf ...
    Pw Psp Psu Pf ...
    Ta2 Sa Ra Pa ...
    Flowsp Flowsu Flowf Flowa]
sz = size(VarMod);

for i = 1:sz(2)
[~, ~, h(i), sig(i), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ktaub([VarMod(:, i), yr'], 0.05, 1)
end 
plot (h)
%% Mid values
[Tm,~, timem] = MonthlyHruWeigthedBasinSum(timeCRHM,Ta, [4,12,14,19], hru_area); 
[Pm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,precip, [4,12,14,19], hru_area); 
[Sm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,snow, [4,12,14,19], hru_area); 
[Rm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,rain, [4,12,14,19], hru_area); 
[Bm, ~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,basinflow/3600, 1, hru_area);

% Seasonal components
for i = 1:length(yr)
   t1 = strcat ('01-Dec-', num2str(yr(i)-1));
   t2 = strcat ('01-Mar-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem ==datetime(t2));
Tw(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sw(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rw(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pw (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flow (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Apr-', num2str(yr(i)));
   t2 = strcat ('01-Jul-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsp(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Ssp(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rsp(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Psp (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowsp (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Jul-', num2str(yr(i)));
   t2 = strcat ('01-Oct-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsu(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Ssu(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rsu(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Psu (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowsu (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Oct-', num2str(yr(i)));
   t2 = strcat ('01-Dec-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tf(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sf(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Rf(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pf (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowf (i, 1) = sum(Bm (a:b)); %m w.e per hru

   t1 = strcat ('01-Oct-', num2str(yr(i)-1));
   t2 = strcat ('01-Oct-', num2str(yr(i)));
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Ta2(i, 1)  = mean(Tm (a:b)); %m w.e per hru
Sa(i, 1)  = sum(Sm (a:b)); %m w.e per hru
Ra(i, 1)  = sum(Rm (a:b)); %m w.e per hru
Pa (i, 1) = sum(Pm (a:b)); %m w.e per hru
Flowa (i, 1) = sum(Bm (a:b)); %m w.e per hru
end 

VarMod = [Tw Tsp Tsu Tf Sw Ssp Ssu Sf Rw Rsp Rsu Rf Pw Psp Psu Pf Ta2 Sa Ra Pa Flowsp Flowsu Flowf Flowa]
sz = size(VarMod);

for i = 1:sz(2)
[H(i), p(i)] = Mann_Kendall (VarMod (:, i), 0.05);
[~, ~, h(i), sig(i), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ktaub([VarMod(:, i), yr'], 0.05, 1)

end 


plot(H)

%% For flow component
load('D:\4_PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWEmelt', 'basinflow', 'basingw', 'firnmelt', 'icemelt', 'infil', 'runoff')
rainfallrunoff = (runoff + infil) -icemelt/24  -firnmelt/24;
rainfallrunoff(rainfallrunoff<0)=0;
rainfallrunoff=rainfallrunoff*1000;
icemelt = icemelt/24 *1000;
firnmelt = firnmelt/24*1000;
swemelt = SWEmelt/24*1000;
timeCRHM = time;
hru_area=hruarea;
hru_elev = hruelev;

[Tm,~, timem] = MonthlyHruWeigthedBasinSum(timeCRHM,rainfallrunoff, 1:37, hru_area); 
[Pm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,swemelt, 1:37, hru_area); 
[Sm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,icemelt, 1:37, hru_area); 
[Rm,~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,firnmelt, 1:37, hru_area); 
[Bm, ~, ~] = MonthlyHruWeigthedBasinSum(timeCRHM,basinflow/3600, 1, hru_area);

varname = {'RainfallRunonff_w' 'RainfallRunonff_s_p' 'RainfallRunonff_s_u' 'RainfallRunonff_f' ...
           'Icemelt_W}' 'Icemelt_s_p' 'Icemelt_s_u' 'Icemelt_f'  ...
           'Firnmelt_w' 'Firnmelt_s_p'  'Firnmelt_s_u' 'Firnmelt_f' ...
           'RainfallRunonff_a' 'Snowmelt_a' 'Icemelt_a' 'Firnmelt_a' ... 
           'Snowmelt_w' 'Snowmelt_s_p' 'Snowmelt_s_u' 'Snowmelt_f' ...
          'Streamflow_s_p' 'Streamflow_s_u' 'Streamflow_f' 'Streamflow_a'};
% Seasonal components
for i = 1:length(yr)
% anuual
   t1 = strcat ('01-Oct-', num2str(yr(i)-1), {' '}, '00:00');
   t2 = strcat ('01-Sep-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Ta2(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Sa(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Ra(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Pa(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flowa (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% Fall (Oct-Dec)
   t1 = strcat ('01-Sep-', num2str(yr(i)-1), {' '}, '00:00');
   t2 = strcat ('01-Nov-', num2str(yr(i)-1), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tf(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Sf(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rf(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Pf(i, 1) = nansum(Pm (a:b)); %m w.e per hru
Flowf (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% Winter (Dec-marc)
    t1 = strcat ('01-Dec-', num2str(yr(i)-1), {' '}, '00:00');
   t2 = strcat ('01-Feb-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem ==datetime(t2));
Tw(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Sw(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rw(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Pw(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flow (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% Spring : Apr-July
   t1 = strcat ('01-Mar-', num2str(yr(i)), {' '}, '00:00');
   t2 = strcat ('01-May-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsp(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Ssp(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rsp(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Psp(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flowsp (i, 1) = nansum(Bm (a:b)); %m w.e per hru

% summer
   t1 = strcat ('01-Jun-', num2str(yr(i)), {' '}, '00:00');
   t2 = strcat ('01-Aug-', num2str(yr(i)), {' '}, '00:00');
a = find(timem==datetime(t1));
b = find(timem==datetime(t2));
Tsu(i, 1)  = nanmean(Tm (a:b)); %m w.e per hru
Ssu(i, 1)  = nansum(Sm (a:b)); %m w.e per hru
Rsu(i, 1)  = nansum(Rm (a:b)); %m w.e per hru
Psu(i, 1)  = nansum(Pm (a:b)); %m w.e per hru
Flowsu (i, 1) = nansum(Bm (a:b)); %m w.e per hru
end 

VarMod= [Tw Tsp Tsu Tf ...
           Sw Ssp Ssu Sf ...
           Rw Rsp Rsu Rf  ...
           Pw Psp  Psu Pf ...
           Ta2 Sa Ra Pa ...
          Flowsp Flowsu Flowf Flowa]
sz = size(VarMod);

% VarMod (10, :) = nan; 

for i = 1:sz(2)
[H(i), p(i)] = Mann_Kendall (VarMod (:, i), 0.05);
[~, ~, h(i), sig(i), ~, ~, sigma(i), sen(i),~, ~, ~, ~, ~, ~, ~, ~] = ktaub([yr', VarMod(:, i)], 0.05, 1)
end
sign_var = find(H);
sign_var=sign_var;

varname_sign = varname (sign_var)
sign_slope = sen(sign_var);
sig_sig =sig(sign_var);
% remove winter rainfall runoff, and firnmelt trend - as they moslty rely
% on initial firn condition and time of firn to snow transition
id = [2,3,6,7]
close all
fig = figure('units','inches','outerposition',[0 0 7 5]);

subplot(2,2, 1);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(2))], 0.05, 1);
ylabel ({'Summer Rainfall Runoff (mm)'})
% legend ('data', 'Sens slope', 'Slope confidence interval (95%)')
ylim([-10 100])

subplot(2,2, 2);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(3))], 0.05, 1);
ylabel ({'Summer Icemelt (mm)'})

subplot(2,2,3);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(6))], 0.05, 1);
ylabel ({'Winter Snowmelt (mm)'})

subplot(2,2, 4);
[~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(7))], 0.05, 1);
ylabel ({'Summer Snowmelt (mm)'})
ylim ([-50 300])

filename = 'SignificantTrends_FlowComp_1989_2020'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))

%%
% all the significant
% close all
% fig = figure('units','inches','outerposition',[0 0 8 6]);
% for i = 1:length(sign_var)
% subplot(2,3, i);
% [~ ,~ ,~ ,~ ,~, ~, ~, ~ ,~ ,senplot, ~, ~] = ktaub([ yr', VarMod(:, sign_var(i))], 0.05, 1);
% ylabel (varname_sign(i))
% end 
% % legend ('data', 'Sens slope', 'Slope confidence interval (95%)')


filename = 'SignificantTrends_Comp_1989_2020'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))
