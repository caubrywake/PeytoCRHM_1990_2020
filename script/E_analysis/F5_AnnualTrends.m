% Figure 5_ Annual Ta, P and streamflow component and F7 scatter plot
% Load CRHM resuls
addpath('D:\PeytoCRHM_1990_2020')
clear all
close all

%% Import CRHM results
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\HRUelev_area.mat')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'SWE', 'hru_rain', 'hru_snow', 'hru_t', 'time')
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210323.mat', 'basinflow', 'basingw')

savedir = 'D:\PeytoCRHM_1990_2020\data_process\chrmoutput\'
figdir = 'D:\PeytoCRHM_1990_2020\fig\analysis\'
 
%% Compile annual values for annual and seasonal hydrometeorological variables
snow = hru_snow;
rain = hru_rain;
precip = snow+rain;
Ta=hru_t;
timeCRHM = time;
hru_area=hruarea;
hru_elev = hruelev;
yr = 1988:2020;
hru = 1:37;
for i = 1:length (yr)
 t1 = strcat ('01-Oct-', num2str(yr(i)-1), {' '}, '01:00');
 t2 = strcat ('30-Sep-', num2str(yr(i)),{' '}, '01:00');
   
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));

x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A(i, 1) = xc;

x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A(i, 2) = xc;

x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
xc = xc-273.15;
A(i, 3) = xc;

A(i, 4)= sum(basinflow(a:b)+basingw(a:b))/(sum(hru_area*10^6));

% Fall
 t1 = strcat ('01-Sep-', num2str(yr(i)-1),{' '}, '01:00');
 t2 = strcat ('30-Nov-', num2str(yr(i)-1),{' '}, '01:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A(i, 5)  = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A(i, 6)  = xc;

x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
A(i, 7) = xc-273.15;

A(i, 8) = sum(basinflow(a:b)+basingw(a:b))/(sum(hru_area*10^6))*1000;

% Winter
 t1 = strcat ('01-Dec-', num2str(yr(i)-1),{' '}, '01:00');
 t2 = strcat ('28-Feb-', num2str(yr(i)),{' '}, '01:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A(i, 9)  = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A(i, 10) = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
A(i, 11) = xc-273.15;
A(i, 12) = sum(basinflow(a:b))/(sum(hru_area*10^6))*1000;

% Spring
 t1 = strcat ('01-Mar-', num2str(yr(i)),{' '}, '01:00');
 t2 = strcat ('30-May-', num2str(yr(i)),{' '}, '01:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A(i, 13)  = xc;
x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A(i, 14)  = xc;
x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
A(i, 15) = xc-273.15;
A(i, 16)  = sum(basinflow(a:b))/(sum(hru_area*10^6))*1000;

% Summer
 t1 = strcat ('01-Jun-', num2str(yr(i)),{' '}, '01:00');
 t2 = strcat ('30-Aug-', num2str(yr(i)),{' '}, '01:00');
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
x = snow;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A(i, 17)  = xc;

x = rain;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A(i, 18)  = xc;

x = Ta+273.15;
xc  = sum(mean(x(a:b,hru)).*(hru_area(hru)/sum(hru_area))); 
A(i, 19) = xc-273.15;
A(i, 20)  = sum(basinflow(a:b))/(sum(hru_area*10^6))*1000;

% peak swe
 t1 = strcat ('01-Oct-', num2str(yr(i)-1),{' '}, '01:00');
 t2 = strcat ('30-Sep-', num2str(yr(i)),{' '}, '01:00');
   
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));
A(i,21) = max(SWE(a:b));
% annual temp at high and low
x = Ta+273.15;
xc  = mean(mean(x(a:b,10))); 
A(i, 22) = xc-273.15;
x = Ta+273.15;
xc  = mean(mean(x(a:b,1))); 
A(i, 23) = xc-273.15;

% snow at low elev
x = snow;
X=  find(hru_elev < 2200);
xa  = sum(x(a:b, X))./1000; xb = xa.*hru_area(X)*10^6 ;xc = sum(xb)/(sum(hru_area(X)*10^6)); % m3 w.e. 
A(i, 24) = xc;

% rain at low elev
x = rain;
xa  = sum(x(a:b, X))./1000; xb = xa.*hru_area(X)*10^6 ;xc = sum(xb)/(sum(hru_area(X)*10^6)); % m3 w.e. 
A(i, 25) = xc;

% total precip
x = precip;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area(hru)*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A(i, 26)  = xc;

% precip at low elev
x = precip;
xa  = sum(x(a:b, X))./1000; xb = xa.*hru_area(X)*10^6 ;xc = sum(xb)/(sum(hru_area(X)*10^6)); % m3 w.e. 
A(i, 29) = xc;

% min temp
x = Ta+273.15;
xc  = mean(min(x(a:b,hru))); 
A(i, 31) = xc-273.15;

% max temp
x = Ta+273.15;
xc  = mean(max(x(a:b,hru))); 
A(i, 32) = xc-273.15;

% snow at hig elev
x = snow;
X=  find(hru_elev > 2700);
xa  = sum(x(a:b, X))./1000; xb = xa.*hru_area(X)*10^6 ;xc = sum(xb)/(sum(hru_area(X)*10^6)); % m3 w.e. 
A(i, 33) = xc;
% precip at hig elev
x = precip;
xa  = sum(x(a:b, X))./1000; xb = xa.*hru_area(X)*10^6 ;xc = sum(xb)/(sum(hru_area(X)*10^6)); % m3 w.e. 
A(i, 34) = xc;
end

A(:, 27) =  A(:, 2)./A(:, 26)*100; %rain/total precip
A(:, 28) =  A(:, 1)./A(:, 26)*100; %snow/ total precip
A(:, 30) = A(:, 24)./A(:, 29)*100; % low elev snow/low elev precip

A(:, 35) = A(:, 33)./A(:, 34)*100; % high elev snow/high elev precip precip

%% make it into a table with headers
A = round(A, 2);

VarCRHMBasin = table(yr', A(:, 1), A(:, 2), A(:, 3), A(:, 4), A(:, 5), A(:, 6),...
   A(:, 7), A(:, 8), A(:, 9), A(:, 10), A(:, 11), A(:, 12),A(:, 13), A(:, 14),...
   A(:, 15), A(:, 16), A(:, 17), A(:, 18), A(:, 19), A(:, 20),A(:, 21), A(:, 22),...
   A(:, 23), A(:, 24), A(:, 25), A(:,26), A(:, 27), A(:, 28),A(:, 29), A(:, 30),A(:, 31), A(:, 32),...
   A(:,33), A(:, 34), A(:, 35));

   lab = { 'year', 'snow','rain','ta','flow', ...
    'snow_fall','rain_fall','ta_fall','flow_fall', ...
    'snow_winter','rain_winter','ta_winter','flow_winter', ...  
    'snow_spring','rain_spring','ta_spring','flow_spring', ...  
    'snow_summer','rain_summer','ta_summer','flow_summer', ...  
    'peakswe','ta_hru10', 'ta_hru1', 'snow_low', ....
    'rain_low', 'total_precip','rainratio', 'snowratio',...
    'low_precip','ratio_lowsnow', 'ta_min','ta_max',...
    'high_snow', 'hig_precip','highsnow_ratio'}
VarCRHMBasin.Properties.VariableNames = lab;
%
writetable(VarCRHMBasin, strcat(savedir, 'VarCRHMBasin.txt'))
save(strcat(savedir, 'VarCRHMBasin.mat'), 'VarCRHMBasin'); 

%% VarCRHMratio_stat 
B (1,:)= max(A);
B (2,:)= min(A);
B (3,:)= nanmean(A);
B (4,:)= nanstd(A);
B = round(B, 2);

varname = {'max';'min';'mean';'std'};
VarCRHMBasin_stat = table(varname,  B(:, 1), B(:, 2), B(:, 3), B(:, 4), B(:, 5), B(:, 6),...
   B(:, 7), B(:, 8), B(:, 9), B(:, 10), B(:, 11), B(:, 12),B(:, 13), B(:, 14),...
   B(:, 15), B(:, 16), B(:, 17), B(:, 18), B(:, 19), B(:, 20),B(:, 21), B(:, 22),...
   B(:, 23), B(:, 24), B(:, 25), B(:,26), B(:, 27), B(:, 28),B(:, 29), B(:, 30),B(:, 31), B(:, 32),...
   B(:,33), B(:, 34), B(:, 35));
   lab = { 'year', 'snow','rain','ta','flow', ...
    'snow_fall','rain_fall','ta_fall','flow_fall', ...
    'snow_winter','rain_winter','ta_winter','flow_winter', ...  
    'snow_spring','rain_spring','ta_spring','flow_spring', ...  
    'snow_summer','rain_summer','ta_summer','flow_summer', ...  
    'peakswe','ta_hru10', 'ta_hru1', 'snow_low', ....
    'rain_low', 'total_precip', 'rainratio','snowratio',...
    'low_precip','ratio_lowsnow', 'ta_min','ta_max',...
    'high_snow', 'hig_precip','highsnow_ratio'}
VarCRHMBasin_stat.Properties.VariableNames = lab;
writetable(VarCRHMBasin_stat, strcat(savedir, 'VarCRHMBasin_stat.txt'))
save(strcat(savedir, 'VarCRHMBasin_stat.mat'), 'VarCRHMBasin_stat'); 


%% Streamflow components
load('D:\PeytoCRHM_1990_2020\data_process\chrmoutput\output_20210314.mat', 'SWEmelt', 'basinflow', 'basingw', 'firnmelt', 'icemelt', 'infil', 'runoff')
rainfallrunoff = (runoff + infil) -icemelt/24  -firnmelt/24;
rainfallrunoff(rainfallrunoff<0)=0;
swemelt = SWEmelt;
yr = 1988:2020
hru = 1:37;
basinflow = basinflow+basingw;
for i = 1:length (yr)
 t1 = strcat ('01-Oct-', num2str(yr(i)-1), {' '}, '01:00');
 t2 = strcat ('30-Sep-', num2str(yr(i)), { ' '}, '01:00');
   
a = find(timeCRHM==datetime(t1));
b = find(timeCRHM==datetime(t2));

x = swemelt/24;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ;xc = sum(xb)/(sum(hru_area*10^6)); % m3 w.e. 
A2(i, 1) = xc;

x = icemelt/24;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A2(i, 2) = xc;

x = firnmelt/24;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A2(i, 3) = xc;

x = rainfallrunoff;
xa  = sum(x(a:b, hru))./1000; xb = xa.*hru_area*10^6 ; xc = sum(xb)/(sum(hru_area*10^6)); 
A2(i, 4) = xc;

A2(i, 5) = sum(basinflow(a:b))/(sum(hru_area*10^6)); %m
end
 
A2(:, 6) = sum(A2(:, 1:4),2)
A2 = round(A2, 2);
FlowComponent_CRHMbasin = table(yr', A2(:, 1), A2(:, 2), A2(:, 3), A2(:, 4), A2(:, 5), A2(:, 6));
   lab = { 'year', 'swemelt', 'icemelt','firnmelt','rainfallrunoff','basinflow', 'basinflow_sumcomp'}
FlowComponent_CRHMbasin.Properties.VariableNames = lab;
writetable(FlowComponent_CRHMbasin, strcat(savedir, 'FlowComponent_CRHMbasin.txt'))
save(strcat(savedir, 'FlowComponent_CRHMbasin.mat'), 'FlowComponent_CRHMbasin'); 

%% Flow component stats
clear B2
B2 (1,:)= max(A2);
B2 (2,:)= min(A2);
B2 (3,:)= nanmean(A2);
B2 (4,:)= nanstd(A2);
B2 = round(B2, 2);

varname = {'max';'min';'mean';'std'};
FlowComponent_CRHMbasin_stat = table(varname, B2(:, 1), B2(:, 2), B2(:, 3), B2(:, 4), B2(:, 5), B2(:, 6));
   lab = { 'year', 'swemelt', 'icemelt','firnmelt','rainfallrunoff','basinflow', 'basinflow_sumcomp'}
FlowComponent_CRHMbasin_stat.Properties.VariableNames = lab;
writetable(FlowComponent_CRHMbasin_stat, strcat(savedir, 'FlowComponent_CRHMbasin_stat.txt'))
save(strcat(savedir, 'FlowComponent_CRHMbasin_stat.mat'), 'FlowComponent_CRHMbasin_stat'); 
%%
% Flow ratios
for i = 1:4;
FlowComponentRatios(:, i) = round(A2(:, i) ./ sum(A2(:,1:4),2).*100)
end
FlowComponentRatio_CRHMbasin = table(yr', FlowComponentRatios(:, 1), FlowComponentRatios(:, 2), FlowComponentRatios(:, 3), FlowComponentRatios(:, 4));
   lab = { 'year', 'swemelt_ratio', 'icemelt_ratio','firnmelt_ratio','rainfallrunoff_ratio'}
FlowComponentRatio_CRHMbasin.Properties.VariableNames = lab;
writetable(FlowComponentRatio_CRHMbasin, strcat(savedir, 'FlowComponentRatio_CRHMbasin.txt'))
save(strcat(savedir, 'FlowComponentRatio_CRHMbasin.mat'), 'FlowComponentRatio_CRHMbasin'); 

%% Flow component stats
B3 (1,:)= max(FlowComponentRatios);
B3 (2,:)= min(FlowComponentRatios);
B3 (3,:)= nanmean(FlowComponentRatios);
B3 (4,:)= nanstd(FlowComponentRatios);
B3 = round(B3, 2);
FlowComponentRatio_CRHMbasin_stat = table(varname, B3(:, 1), B3(:, 2), B3(:, 3), B3(:, 4));
   lab = { 'year', 'swemelt_ratio', 'icemelt_ratio','firnmelt_ratio','rainfallrunoff_ratio'}
FlowComponentRatio_CRHMbasin_stat.Properties.VariableNames = lab;
writetable(FlowComponentRatio_CRHMbasin_stat, strcat(savedir, 'FlowComponentRatio_CRHMbasin_stat.txt'))
save(strcat(savedir, 'FlowComponentRatio_CRHMbasin_stat.mat'), 'FlowComponentRatio_CRHMbasin_stat'); 

%% Figure with the summer and annual temperature, rain and snow and basinflow
% COLOR PALETTE
cavy = [102 170 255]/255%pale blue
csnow= [0 102 204]/255%  mid blue
crain = [0 51 102]/255% dark blue

crunoff = [0 153 153]/255 % light green
cdrift = [0 102 102]/255% dark green
cET = [0 204 204]/255 % black

csnowm = [200 200 200]/255% light grey
cicem =  [140 140 140]/255% mid grey 
cfirnm = [90 90 90]/255% dark grey

% MAIN FIGURE
fig = figure('units','inches','outerposition',[0 0 8 5.5]);
yr = 1988:2020
subplot (3,1,1) %air T
plot (yr, A(:,3) , '-xk'); hold on
ylabel ({'Temperature'; '(^{\circ}C)'})
ylim ([-5.5 -2])
yticks ([-5:1:-2])
xlim ([1987.5 2020.5])
text(1987.5, -1.7, 'a)')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot   (3,1,2) % precip (snow and rain)
h = bar(yr, [A(:,1), A(:,2)],'stack')%, 'BarWidth', 1.1)
set(h(1), 'FaceColor', csnow)
set(h(2), 'FaceColor', crain)
xlim ([1987.5 2020.5])
ylim ([0 3.2])
grid on
box on
ylabel ({'Precipitation';'(m)'})
text(1987.5, 3.55, 'b)')
h = legend ('Snow', 'Rain', 'Orientation', 'Horizontal', 'Location', 'Southwest')
pos = get(h,'Position');
set(h,'Position',[pos(1)-0.01 pos(2)-0.01 pos(3) pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Streamflow Component
subplot   (3,1,3) % precip (snow and rain)
h = bar(yr, [A2(:,1) A2(:,2) A2(:,3) A2(:,4)],'stack')%, 'BarWidth', 1.1)
set(h(1), 'FaceColor', csnowm)
set(h(2), 'FaceColor', cicem)
set(h(3), 'FaceColor', cfirnm)
set(h(4), 'FaceColor',crunoff)
xlim ([1987.5 2020.5])
ylim ([0 3.5])
grid on
box on
ylabel ({'Streamflow';'Component';'(m. w.e.)'})
text(1987.5, 3.95, 'c)')
h2 = legend ('Snowmelt', 'Icemelt', 'Firnmelt', 'Rainfall Runoff', 'Orientation', 'Horizontal', 'Location', 'Southwest')
pos = get(h2,'Position');
set(h2,'Position',[pos(1)-0.01 pos(2)-0.01 pos(3) pos(4)]);
tightfig(fig)
% 
%
filename = 'AnnualValues_1987_2020'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))




